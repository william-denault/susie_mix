"""Local fixtures only: Slurm calls are mocked; no jobs or real results touched."""
import contextlib
import csv
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest

PROJECT = Path(__file__).resolve().parents[3]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


recovery = module("recovery", PROJECT / "script/sim/recover_failed_checkpoints.py")
shared = module("shared", PROJECT / "script/scan_tissue_attempt/slide_shared_scheduler.py")


class Scheduler:
    paths = staticmethod(shared.paths)
    read = staticmethod(shared.read)
    digest = staticmethod(shared.digest)
    check_config = staticmethod(shared.check_config)
    atomic = staticmethod(shared.atomic)

    def __init__(self):
        self.rows = [dict(id="9000", array="9000", name="slide_parallel_control",
                          state="PENDING", reason="BeginTime")]
        self.accounting = "\n".join("{}_{}|FAILED|1:0".format(recovery.ARRAY, t) for t in recovery.TASKS)
        self.calls = []

    @contextlib.contextmanager
    def lock(self, path):
        yield

    def queue(self):
        return [dict(row) for row in self.rows]

    def run(self, args):
        self.calls.append(list(args))
        if args[0] == "sacct":
            return self.accounting
        if args[:2] == ["scontrol", "hold"]:
            self.rows[0]["reason"] = "JobHeldUser"
            return ""
        if args[:2] == ["scontrol", "release"]:
            self.rows[0]["reason"] = "BeginTime"
            return ""
        raise AssertionError("Unexpected external command: " + repr(args))


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        base = (PROJECT / "tmp/n250_simulation_validation").resolve()
        base.mkdir(parents=True, exist_ok=True)
        self.directory = tempfile.TemporaryDirectory(prefix="checkpoint_recovery_", dir=str(base))
        self.project = Path(self.directory.name).resolve()
        assert self.project.parent == base  # Bound automatic recursive fixture cleanup.
        self.addCleanup(self.directory.cleanup)
        self.scheduler = Scheduler()
        self.root = shared.paths(self.project)
        self.receipt = self.root / ("checkpoint_recovery_" + recovery.ARRAY + ".json")
        for name in shared.SIM_SOURCES:
            f = self.project / name
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_text("frozen source\n")
        for batch in (2, 3):
            f = self.project / "job/run_simulation_slide_batch_{}".format(batch)
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_text(shared.bridge(batch))
        self.config = {"schema": 1, "project": str(self.project), "jobs": list(range(301, 901)),
                       "source_hashes": {name: shared.digest(self.project / name) for name in shared.SIM_SOURCES},
                       "bridges": [2, 3]}
        shared.atomic(self.root / "enabled.json", self.config)
        rows = []
        self.checkpoints = []
        for job in recovery.JOBS:
            f = self.project / "script/sim/jobs_slide/sim_job_{}.R".format(job)
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_text("driver {}\n".format(job))
            rel = "simulation results/n250_r300_v1/chunks/condition_{}.RData".format(job)
            rows.append(dict(job_id=job, n=250, reps_per_chunk=300, schema_version=4, output_file=rel))
            checkpoint = self.project / rel
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            checkpoint.write_bytes(("old checkpoint {}".format(job)).encode())
            self.checkpoints.append(checkpoint)
            log = self.project / "job/logs/susie_slide_sim_{}_{}.err".format(recovery.ARRAY, job - 1)
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text(recovery.MESSAGE + "\nExecution halted\n")
        with (self.project / "script/sim/jobs_slide/manifest.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        self.unrelated = self.checkpoints[0].parent / "unrelated.RData"
        self.unrelated.write_bytes(b"keep successful data")

    def run_recovery(self, apply=True, **kwargs):
        with contextlib.redirect_stdout(io.StringIO()):
            recovery.recover(self.project, self.scheduler, apply=apply, **kwargs)

    def test_dry_run_preserves_everything(self):
        self.run_recovery(apply=False)
        self.assertFalse(self.receipt.exists())
        self.assertTrue(all(p.exists() for p in self.checkpoints))
        self.assertEqual(shared.read(self.root / "enabled.json"), self.config)
        self.assertFalse(any(c[0] == "scontrol" for c in self.scheduler.calls))

    def test_archive_register_and_schedule_exactly_twenty(self):
        temporary = Path(str(self.checkpoints[0]) + ".tmp")
        temporary.write_bytes(b"interrupted save")
        self.run_recovery()
        config = shared.read(self.root / "enabled.json")
        shared.check_config(self.project, config)
        pending = shared.sim_pending(self.project, config, {"arrays": []}, self.scheduler.queue())
        self.assertEqual(pending[:20], list(recovery.JOBS))
        self.assertEqual(len(config["jobs"]), 620)
        self.assertEqual(self.unrelated.read_bytes(), b"keep successful data")
        for p in self.checkpoints:
            self.assertFalse(p.exists())
            backup = self.project / "simulation results/n250_r300_v1/checkpoint_backups/array_49451873" / p.name
            self.assertTrue(backup.read_bytes().startswith(b"old checkpoint"))
        self.assertFalse(temporary.exists())
        self.assertEqual(self.scheduler.rows[0]["reason"], "BeginTime")
        self.assertFalse(any(c[0] == "sbatch" for c in self.scheduler.calls))

    def test_repeated_command_preserves_new_checkpoint_and_completed_receipt(self):
        self.run_recovery()
        self.checkpoints[0].write_bytes(b"new successful results")
        shared.atomic(self.root / "done/{}.json".format(recovery.JOBS[0]), {"job": recovery.JOBS[0]})
        self.run_recovery()
        self.assertEqual(self.checkpoints[0].read_bytes(), b"new successful results")

    def test_interrupted_archive_can_resume(self):
        calls = []
        def interrupted(source, target):
            if calls:
                raise RuntimeError("simulated interruption")
            calls.append(source)
            recovery.os.rename(source, target)
        with self.assertRaisesRegex(RuntimeError, "simulated interruption"):
            self.run_recovery(rename=interrupted)
        self.assertEqual(shared.read(self.root / "enabled.json"), self.config)
        self.assertEqual(self.scheduler.rows[0]["reason"], "BeginTime")
        self.run_recovery()
        self.assertEqual(shared.read(self.receipt)["phase"], "registered")

    def test_interruption_after_registration_does_not_move_new_results(self):
        self.run_recovery()
        receipt = shared.read(self.receipt)
        receipt["phase"] = "prepared"
        receipt["released"] = []
        shared.atomic(self.receipt, receipt)
        self.scheduler.rows[0]["reason"] = "JobHeldUser"
        self.checkpoints[0].write_bytes(b"new run")
        self.run_recovery()
        self.assertEqual(self.checkpoints[0].read_bytes(), b"new run")
        self.assertEqual(self.scheduler.rows[0]["reason"], "BeginTime")

    def test_other_failure_causes_are_rejected_before_moves(self):
        log = self.project / "job/logs/susie_slide_sim_49451873_16.err"
        log.write_text("OUT_OF_MEMORY\n")
        with self.assertRaisesRegex(RuntimeError, "mismatch not confirmed"):
            self.run_recovery()
        self.assertTrue(all(p.exists() for p in self.checkpoints))
        self.assertFalse(self.receipt.exists())

    def test_completed_task_is_not_restarted(self):
        self.scheduler.accounting = self.scheduler.accounting.replace("49451873_16|FAILED|1:0", "49451873_16|COMPLETED|0:0")
        with self.assertRaisesRegex(RuntimeError, "not confirmed FAILED"):
            self.run_recovery()
        self.assertTrue(all(p.exists() for p in self.checkpoints))

    def test_running_target_is_not_moved(self):
        self.scheduler.rows.append(dict(id="49451873_16", array="49451873", name="susie_slide_sim", state="RUNNING", reason="None"))
        with self.assertRaisesRegex(RuntimeError, "still active"):
            self.run_recovery()
        self.assertTrue(all(p.exists() for p in self.checkpoints))

    def test_user_hold_is_preserved(self):
        self.scheduler.rows[0]["reason"] = "JobHeldUser"
        with self.assertRaisesRegex(RuntimeError, "held by someone else"):
            self.run_recovery()
        self.assertFalse(any(c[0] == "scontrol" for c in self.scheduler.calls))
        self.assertTrue(all(p.exists() for p in self.checkpoints))

    def test_source_drift_is_not_bypassed(self):
        (self.project / shared.SIM_SOURCES[0]).write_text("changed math")
        with self.assertRaisesRegex(RuntimeError, "Frozen simulation source changed"):
            self.run_recovery()
        self.assertFalse(self.receipt.exists())


if __name__ == "__main__":
    unittest.main()
