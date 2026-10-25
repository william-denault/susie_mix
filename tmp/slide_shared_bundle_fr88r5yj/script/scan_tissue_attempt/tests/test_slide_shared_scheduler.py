"""Mock Slurm integration tests. No real submissions or production data reads."""
import contextlib
import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location("shared", REPO / "script/scan_tissue_attempt/slide_shared_scheduler.py")
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)


@contextlib.contextmanager
def no_lock(path):
    yield


class SchedulerTest(unittest.TestCase):
    def setUp(self):
        (REPO / "tmp").mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(prefix="shared_test_", dir=str(REPO / "tmp"))
        self.project = Path(self.temp.name).resolve() / "project with spaces"
        self.root = s.paths(self.project)
        self.rows = [self.row("10", "susie_slide_launch"), self.row("20", "slide_parallel_control")]
        self.calls = []
        self.next_id = 1000
        self.em_status = "batch"
        self.em_tasks = None
        self.fail_submit = False
        self.ambiguous = False
        self.write("script/sim/run_job.R", "# unchanged runner\n")
        for source in s.SIM_SOURCES[1:]:
            self.write(source, "# unchanged shared simulation source\n")
        self.write("script/sim/jobs_slide/submission_batches.csv",
                   "batch,array_start,array_end,offset,first_job,last_job,launcher\n" +
                   "".join("{b},0,299,{o},{f},{l},job/run_simulation_slide_batch_{b}\n".format(
                       b=b, o=(b-1)*300, f=(b-1)*300+1, l=b*300) for b in range(1, 4)))
        for b in range(1, 4):
            self.write("job/run_simulation_slide_batch_" + str(b), "# original batch {}\n#SBATCH --array=0-299\n".format(b))
        for job in range(1, 901):
            self.write("script/sim/jobs_slide/sim_job_{}.R".format(job), "# original condition {}\n".format(job))
        self.write("simulation results/n250_r300_v1/launcher/last_batch.txt", "1\n")
        self.write("simulation results/n250_r300_v1/launcher/last_array_job_id.txt", "9\n")
        self.write("simulation results/n250_r300_v1/launcher/last_continuation_job_id.txt", "10\n")
        self.write("simulation results/n250_r300_v1/launcher/submissions.tsv", "1\tarray\t9\n2\tcontinuation\t10\n")
        self.write("results_slide_prior_em/iteration_000/results/sentinel.rds", "DO NOT MODIFY\n")
        self.patches = [patch.object(s, "run", self.mock_run), patch.object(s, "lock", no_lock),
                        patch.dict(os.environ, {"USER": "fixture", "SLURM_JOB_ID": "20",
                                               "SUSIE_SLIDE_WORKERS": "298", "SUSIE_SLIDE_POLL_MINUTES": "5",
                                               "SUSIE_MIX_SIM_INPUT": "gtex", "SUSIE_MIX_SIM_GENE": "",
                                               "SUSIE_SLIDE_WAIT_FOR_SIMULATIONS": "1"})]
        for p in self.patches:
            p.start()

    def tearDown(self):
        for p in reversed(self.patches):
            p.stop()
        self.temp.cleanup()

    def write(self, name, text):
        path = self.project / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    @staticmethod
    def row(job, name, state="PENDING", array=None):
        return {"id": job, "array": array or job.split("_")[0], "name": name,
                "state": state, "reason": "Dependency"}

    def mock_run(self, args):
        args = [str(a) for a in args]
        self.calls.append(args)
        if args[0] == "squeue":
            return "\n".join("{id}|{array}|{name}|{state}|{reason}".format(**r) for r in self.rows)
        if args[0] == "scontrol":
            for row in self.rows:
                if row["id"] == args[2]:
                    row["reason"] = "JobHeldUser" if args[1] == "hold" else "Dependency"
            return ""
        if args[0] == "Rscript":
            cap = int(args[-2])
            self.capacity = cap
            if self.em_status != "batch" or cap == 0:
                return ("waiting" if self.em_status == "batch" else self.em_status) + "\n1"
            batch = self.root.parent / "iteration_001/recovery/parallel/batches/fixture"
            batch.mkdir(parents=True, exist_ok=True)
            tasks = min(cap, self.em_tasks) if self.em_tasks is not None else cap
            return "batch\n1\n{}\n{}\nfit".format(batch, tasks)
        if args[0] == "sbatch":
            if self.fail_submit:
                raise subprocess.CalledProcessError(1, args)
            if self.ambiguous:
                return "response lost"
            self.next_id += 1
            return str(self.next_id) + ";test_cluster"
        raise AssertionError(args)

    def activate(self, simulations=0, legacy=True):
        self.rows += [self.row("9_" + str(i), "susie_slide_sim") for i in range(simulations)]
        s.enable(self.project)
        for row in self.rows:
            if row["id"] == "20":
                row["state"] = "RUNNING"
        if not legacy:
            self.rows = [r for r in self.rows if r["name"] != "susie_slide_launch"]
        self.calls.clear()

    def submissions(self):
        return [a for a in self.calls if a[0] == "sbatch"]

    def additions(self):
        return sum(int(a.split("-")[-1]) for call in self.submissions() for a in call if a.startswith("--array="))

    def test_handover_preserves_existing_and_reuses_chain(self):
        old = (self.project / "job/run_simulation_slide_batch_1").read_bytes()
        self.activate(simulations=290)
        config = s.read(self.root / "enabled.json")
        self.assertEqual(config["jobs"], list(range(301, 901)))
        self.assertEqual(old, (self.project / "job/run_simulation_slide_batch_1").read_bytes())
        self.assertIn("# original", (self.root / "original_launchers/run_simulation_slide_batch_2").read_text())
        self.assertNotIn("#SBATCH --array", (self.project / "job/run_simulation_slide_batch_2").read_text())
        self.assertEqual((self.project / "results_slide_prior_em/iteration_000/results/sentinel.rds").read_text(), "DO NOT MODIFY\n")
        self.assertFalse(self.submissions())  # No second controller on activation.
        self.assertTrue(all(r["reason"] == "Dependency" for r in self.rows))
        self.calls.clear()
        s.enable(self.project)  # Idempotent, no duplicate requests/submissions.
        self.assertFalse(self.submissions())

    def test_full_budget_balances_149_each_and_five_minutes(self):
        self.activate(legacy=False)
        s.control(self.project)
        self.assertEqual(self.capacity, 149)
        self.assertEqual(self.additions(), 298)
        self.assertEqual(len(s.read(self.root / "state.json")["arrays"][0]["jobs"]), 149)
        self.assertIn("--begin=now+5minutes", self.submissions()[-1])
        self.assertIn("--dependency=afterany:20", self.submissions()[-1])

    def test_existing_simulations_free_slots_go_to_em_no_cancel(self):
        self.activate(simulations=250)
        s.control(self.project)
        self.assertEqual(self.capacity, 45)  # self + legacy + next + two legacy reserve
        self.assertEqual(self.additions(), 45)
        self.assertEqual(len(self.submissions()), 2)  # EM + continuation
        self.assertFalse(any(a[0] == "scancel" for a in self.calls))

    def test_one_freed_slot_refilled_and_other_user_jobs_counted(self):
        self.activate(simulations=148, legacy=False)
        self.rows += [self.row("30_" + str(i), "slide_parallel_fit") for i in range(147)]
        self.rows += [self.row("40_" + str(i), "unrelated") for i in range(2)]
        s.control(self.project)
        self.assertEqual(self.capacity, 1)
        self.assertEqual(self.additions(), 1)
        self.assertLessEqual(len(self.rows) + self.additions() + 1, 300)

    def test_em_completion_lends_all_slots_to_simulations_and_keeps_chain(self):
        self.activate(legacy=False)
        self.em_status = "complete"
        s.control(self.project)
        self.assertEqual(self.additions(), 298)
        self.assertEqual(len(self.submissions()), 2)

    def test_waiting_em_lends_unused_capacity(self):
        self.activate(legacy=False)
        self.em_tasks = 3
        s.control(self.project)
        self.assertEqual(len(s.read(self.root / "state.json")["arrays"][0]["jobs"]), 295)
        self.assertEqual(self.additions(), 298)

    def test_simulations_finished_em_can_use_298(self):
        self.activate(legacy=False)
        config = s.read(self.root / "enabled.json")
        config["jobs"] = []
        s.atomic(self.root / "enabled.json", config)
        s.control(self.project)
        self.assertEqual(self.capacity, 298)
        self.assertEqual(self.additions(), 298)

    def test_duplicate_controller_and_unknown_sim_worker_stop(self):
        self.activate(legacy=False)
        self.rows.append(self.row("21", "slide_parallel_control"))
        with self.assertRaisesRegex(RuntimeError, "Another EM"):
            s.control(self.project)
        self.rows.pop()
        self.rows.append(self.row("99_1", s.SIM_WORKER))
        with self.assertRaisesRegex(RuntimeError, "Unregistered"):
            s.control(self.project)
        self.assertFalse(self.submissions())

    def test_ambiguous_submission_blocks_relaunch(self):
        self.activate(legacy=False)
        self.ambiguous = True
        with self.assertRaisesRegex(RuntimeError, "Ambiguous"):
            s.control(self.project)
        self.calls.clear()
        self.ambiguous = False
        with self.assertRaisesRegex(RuntimeError, "Unresolved"):
            s.control(self.project)
        self.assertFalse(self.submissions())

    def test_migration_running_or_preheld_launchers_untouched(self):
        self.rows[0]["state"] = "RUNNING"
        with self.assertRaisesRegex(RuntimeError, "running"):
            s.enable(self.project)
        self.rows[0]["state"] = "PENDING"
        self.rows[0]["reason"] = "JobHeldUser"
        with self.assertRaisesRegex(RuntimeError, "already held"):
            s.enable(self.project)
        self.assertFalse(any(a[0] == "scontrol" for a in self.calls))

    def test_partial_handover_resumes_and_releases_own_holds(self):
        original = s.atomic
        def fail_enable(path, *args, **kwargs):
            if Path(path).name == "enabled.json":
                raise OSError("interrupted activation")
            return original(path, *args, **kwargs)
        with patch.object(s, "atomic", fail_enable):
            with self.assertRaisesRegex(OSError, "interrupted"):
                s.enable(self.project)
        self.assertTrue(all(r["reason"] == "JobHeldUser" for r in self.rows))
        s.enable(self.project)
        self.assertTrue(all(r["reason"] == "Dependency" for r in self.rows))
        self.assertEqual(s.read(self.root / "enabled.json")["jobs"], list(range(301, 901)))

    def test_worker_resume_receipt_and_failure_budget(self):
        self.activate(legacy=False)
        config = s.read(self.root / "enabled.json")
        config["jobs"] = [301]
        s.atomic(self.root / "enabled.json", config)
        batch = self.root / "batches/000001"
        s.atomic(batch / "assignment.json", {"jobs": [301]})
        state = {"arrays": [{"array": "50", "batch": str(batch), "jobs": [301]}]}
        self.assertEqual(s.sim_pending(self.project, config, state, [self.row("50_1", s.SIM_WORKER)]), [])
        self.assertEqual(s.sim_pending(self.project, config, state, []), [301])
        with patch.object(s.subprocess, "call", return_value=0) as call:
            s.worker(self.project, batch, 1)
            self.assertEqual(call.call_args[0][0][-1], str(self.project / "script/sim/jobs_slide/sim_job_301.R"))
        self.assertEqual(s.sim_pending(self.project, config, state, []), [])
        (self.root / "done/301.json").unlink()
        state["arrays"] *= 3
        with self.assertRaisesRegex(RuntimeError, "three attempts"):
            s.sim_pending(self.project, config, state, [])

    def test_changed_frozen_source_or_bridge_blocks_refill(self):
        self.activate(legacy=False)
        self.write("script/sim/jobs_slide/sim_job_301.R", "changed")
        with self.assertRaisesRegex(RuntimeError, "Frozen simulation source"):
            s.control(self.project)
        self.assertFalse(self.submissions())

    def test_later_batch_already_submitted_is_preserved(self):
        self.write("simulation results/n250_r300_v1/launcher/last_batch.txt", "2\n")
        self.write("simulation results/n250_r300_v1/launcher/last_array_job_id.txt", "11\n")
        self.write("simulation results/n250_r300_v1/launcher/submissions.tsv", "1\tarray\t9\n2\tarray\t11\n3\tcontinuation\t10\n")
        old = (self.project / "job/run_simulation_slide_batch_2").read_bytes()
        self.rows.append(self.row("11_1", "susie_slide_sim"))
        self.activate()
        self.assertEqual(s.read(self.root / "enabled.json")["jobs"], list(range(601, 901)))
        self.assertEqual(old, (self.project / "job/run_simulation_slide_batch_2").read_bytes())

    def test_failed_submission_and_failed_worker_preserve_retry_state(self):
        self.activate(legacy=False)
        self.fail_submit = True
        with self.assertRaises(subprocess.CalledProcessError):
            s.control(self.project)
        self.assertFalse((self.root.parent / "parallel_pending_submission.txt").exists())
        batch = self.root / "batches/000001"
        s.atomic(batch / "assignment.json", {"jobs": [301]})
        with patch.object(s.subprocess, "call", return_value=1):
            with self.assertRaisesRegex(RuntimeError, "exited 1"):
                s.worker(self.project, batch, 1)
        self.assertFalse((self.root / "done/301.json").exists())
        self.assertEqual(s.read(batch / "exit_1.json")["exit_code"], 1)

    def test_inherited_overfull_queue_has_no_new_worker_submissions(self):
        self.activate(simulations=300)
        s.control(self.project)
        self.assertEqual(self.additions(), 0)
        self.assertEqual(len(self.submissions()), 1)  # keep only the small monitor alive


if __name__ == "__main__":
    unittest.main(verbosity=2)
