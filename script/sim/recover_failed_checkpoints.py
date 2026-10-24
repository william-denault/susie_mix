"""Preserve incompatible checkpoints and adopt the 20 failed legacy conditions.

Python >= 3.6. Uses the existing shared scheduler's locks and frozen-source checks.
Default is read-only; --apply archives only these checkpoints and registers retries.
No sbatch, cancellation, package changes, source changes or compatibility bypass.
"""
import argparse
import csv
import importlib.util
import os
from pathlib import Path
import sys

ARRAY = "49451873"
TASKS = tuple(list(range(16, 20)) + list(range(26, 30)) +
              list(range(42, 48)) + list(range(64, 68)) + [70, 71])
JOBS = tuple(t + 1 for t in TASKS)
MESSAGE = "Checkpoint design, input path or package versions differ; use a new output directory."


def load_scheduler(project):
    path = project / "script/scan_tissue_attempt/slide_shared_scheduler.py"
    spec = importlib.util.spec_from_file_location("checkpoint_recovery_scheduler", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_failures(project, scheduler):
    accounting = scheduler.run(["sacct", "-n", "-P", "-X", "-j", ARRAY,
                                "--format=JobID%64,State%32,ExitCode"])
    states = {}
    for line in accounting.splitlines():
        parts = [p.strip() for p in line.split("|")]
        if len(parts) != 3:
            raise RuntimeError("Unexpected sacct record: " + line)
        if parts[0] in states:
            raise RuntimeError("Multiple accounting records for " + parts[0])
        states[parts[0]] = parts[1:]
    for task in TASKS:
        task_id = ARRAY + "_" + str(task)
        if states.get(task_id) != ["FAILED", "1:0"]:
            raise RuntimeError(task_id + " is not confirmed FAILED with exit 1:0; nothing archived")
        log = project / "job/logs" / ("susie_slide_sim_" + task_id + ".err")
        tail = log.read_text(encoding="utf-8", errors="replace")[-16000:]
        if MESSAGE not in tail or "Execution halted" not in tail:
            raise RuntimeError("Inspect " + str(log) + ": checkpoint mismatch not confirmed; nothing archived")


def planned_files(project):
    with (project / "script/sim/jobs_slide/manifest.csv").open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    by_job = {}
    for row in rows:
        job = int(row["job_id"])
        if job in by_job:
            raise RuntimeError("Duplicate manifest job " + str(job))
        by_job[job] = row
    chunk_dir = (project / "simulation results/n250_r300_v1/chunks").resolve()
    backup_dir = project / "simulation results/n250_r300_v1/checkpoint_backups" / ("array_" + ARRAY)
    chunk_dir.relative_to(project)
    backup_dir.resolve().relative_to(project)
    files = []
    for job in JOBS:
        row = by_job[job]
        if (int(row["n"]) != 250 or int(row["reps_per_chunk"]) != 300 or
                int(row["schema_version"]) != 4):
            raise RuntimeError("Unexpected simulation design for condition " + str(job))
        path = project / row["output_file"]
        if path.is_symlink() or path.resolve().parent != chunk_dir or path.suffix != ".RData":
            raise RuntimeError("Checkpoint is outside the intended chunk directory: " + str(path))
        for source in (path, Path(str(path) + ".tmp")):
            if source.is_symlink() or source.resolve().parent != chunk_dir:
                raise RuntimeError("Unexpected checkpoint path: " + str(source))
            target = backup_dir / source.name
            if target.exists():
                raise RuntimeError("Unregistered backup already exists: " + str(target))
            if source.exists() and not source.is_file():
                raise RuntimeError("Checkpoint is not a regular file: " + str(source))
            files.append({"job": job, "source": str(source.relative_to(project)),
                          "backup": str(target.relative_to(project)), "existed": source.exists()})
    return files


def archive_files(project, records, rename=os.rename):
    # A saved transaction makes interruption between individual moves recoverable.
    chunks = (project / "simulation results/n250_r300_v1/chunks").resolve()
    backups = (project / "simulation results/n250_r300_v1/checkpoint_backups" / ("array_" + ARRAY)).resolve()
    chunks.relative_to(project)
    backups.relative_to(project)
    for record in records:
        source, target = project / record["source"], project / record["backup"]
        if (source.is_symlink() or target.is_symlink() or source.resolve().parent != chunks or
                target.resolve().parent != backups):
            raise RuntimeError("Unsafe path in checkpoint recovery transaction")
        if record["existed"]:
            if source.exists() and target.exists():
                raise RuntimeError("Both original and backup exist; inspect before retrying: " + str(source))
            if not source.exists() and not target.exists():
                raise RuntimeError("Original and backup are missing: " + str(source))
            if source.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                rename(str(source), str(target))
        elif source.exists() or target.exists():
            raise RuntimeError("A checkpoint appeared after planning; inspect before retrying: " + str(source))


def release_controllers(scheduler, transaction, receipt):
    live = {r["id"]: r for r in scheduler.queue()}
    for job in transaction.get("hold_ids", []):
        if job in transaction.get("released", []):
            continue
        if job in live and live[job]["state"] == "PENDING" and live[job]["reason"].startswith("JobHeld"):
            scheduler.run(["scontrol", "release", job])
        transaction.setdefault("released", []).append(job)
        scheduler.atomic(receipt, transaction)


def recover(project, scheduler, apply=False, rename=os.rename):
    project = project.resolve()
    root = scheduler.paths(project)
    receipt = root / ("checkpoint_recovery_" + ARRAY + ".json")
    with scheduler.lock(root / ".migration.lock"):
        if (root.parent / "parallel_pending_submission.txt").exists():
            raise RuntimeError("A scheduler submission needs reconciliation; no checkpoints changed")
        config = scheduler.read(root / "enabled.json")
        scheduler.check_config(project, config)
        transaction = scheduler.read(receipt) if receipt.exists() else None
        included = set(JOBS).intersection(config["jobs"])
        if transaction is not None:
            if (transaction.get("array") != ARRAY or transaction.get("jobs") != list(JOBS) or
                    transaction.get("project") != str(project)):
                raise RuntimeError("Recovery receipt belongs to a different operation")
            for name, expected in transaction["source_hashes"].items():
                if scheduler.digest(project / name) != expected:
                    raise RuntimeError("Recovery driver changed: " + name)
            if included == set(JOBS):
                # Registration may have succeeded immediately before interruption.
                # Never archive newly generated results on a repeated invocation.
                if any(config["source_hashes"].get(name) != value
                       for name, value in transaction["source_hashes"].items()):
                    raise RuntimeError("Registered retry hashes do not match the receipt")
                if apply and transaction["phase"] != "registered":
                    transaction["phase"] = "registered"
                    scheduler.atomic(receipt, transaction)
                if apply:
                    release_controllers(scheduler, transaction, receipt)
                print("All 20 retries are already registered. Existing/new results left untouched.")
                return
            if transaction["phase"] != "prepared":
                raise RuntimeError("Recovery registration was removed; inspect the receipt")
        if included:
            raise RuntimeError("One or more target conditions are already owned by the scheduler")
        rows = scheduler.queue()
        active_ids = {r["id"] for r in rows}
        if any(ARRAY + "_" + str(t) in active_ids for t in TASKS):
            raise RuntimeError("A target legacy job is still active; no checkpoints changed")
        controllers = [r for r in rows if r["name"] == "slide_parallel_control"]
        if len(controllers) != 1:
            raise RuntimeError("Expected one active shared controller; no checkpoints changed")
        state = scheduler.read(root / "state.json") if (root / "state.json").exists() else {"arrays": []}
        assigned = {job for array in state["arrays"] for job in array["jobs"]}
        if assigned.intersection(JOBS) or any((root / "done/{}.json".format(j)).exists() for j in JOBS):
            raise RuntimeError("Target conditions already have shared assignments or receipts")
        validate_failures(project, scheduler)
        if transaction is None:
            hashes = {"script/sim/jobs_slide/sim_job_{}.R".format(j):
                      scheduler.digest(project / "script/sim/jobs_slide/sim_job_{}.R".format(j)) for j in JOBS}
            transaction = {"array": ARRAY, "jobs": list(JOBS), "project": str(project),
                           "phase": "prepared", "files": planned_files(project),
                           "source_hashes": hashes, "previous_config": config,
                           "hold_ids": [], "released": []}
        print("Verified all 20 original tasks failed with the checkpoint-mismatch message.")
        print("Conditions: " + ", ".join(str(j) for j in JOBS))
        print("Old checkpoints will be preserved in simulation results/n250_r300_v1/checkpoint_backups/array_" + ARRAY)
        if not apply:
            print("Read-only plan. Add --apply to archive these checkpoints and register fresh runs.")
            return
        if not receipt.exists():
            scheduler.atomic(receipt, transaction)  # Durable intent before the first move.
        # Hold only the pending continuation while editing its job inventory.
        # Record ownership before holding so forced interruption is recoverable.
        controllers = [r for r in scheduler.queue() if r["name"] == "slide_parallel_control"]
        if len(controllers) != 1 or controllers[0]["state"] != "PENDING":
            raise RuntimeError("Controller is currently running; retry this command in five minutes")
        controller = controllers[0]
        owns_hold = controller["id"] in transaction["hold_ids"] and controller["id"] not in transaction["released"]
        if controller["reason"].startswith("JobHeld") and not owns_hold:
            raise RuntimeError("Controller was already held by someone else; no checkpoints changed")
        if controller["id"] not in transaction["hold_ids"]:
            transaction["hold_ids"].append(controller["id"])
        if controller["id"] in transaction["released"]:
            transaction["released"].remove(controller["id"])
        scheduler.atomic(receipt, transaction)
        try:
            scheduler.run(["scontrol", "hold", controller["id"]])
            held = {r["id"]: r for r in scheduler.queue()}.get(controller["id"])
            if held is None or held["state"] != "PENDING" or not held["reason"].startswith("JobHeld"):
                raise RuntimeError("Controller started before it could be held; no checkpoints changed")
            with scheduler.lock(root.parent / ".prepare.lock"):
                current_config = scheduler.read(root / "enabled.json")
                if current_config != config:
                    raise RuntimeError("Scheduler inventory changed; retry after inspecting the recovery receipt")
                if (root.parent / "parallel_pending_submission.txt").exists():
                    raise RuntimeError("Scheduler submission needs reconciliation; no checkpoints changed")
                archive_files(project, transaction["files"], rename=rename)
                updated = dict(config)
                updated["jobs"] = sorted(set(config["jobs"]).union(JOBS))
                updated["source_hashes"] = dict(config["source_hashes"])
                updated["source_hashes"].update(transaction["source_hashes"])
                scheduler.atomic(root / "enabled.json", updated)
                transaction["phase"] = "registered"
                scheduler.atomic(receipt, transaction)
        finally:
            release_controllers(scheduler, transaction, receipt)
        print("REGISTERED: the existing controller will schedule all 20 retries as slots become available.")
        print("No new controller or array submitted. Frozen simulation sources and successful conditions unchanged.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true", help="Preserve conflicting checkpoints and register the 20 retries")
    args = parser.parse_args()
    project = Path(__file__).resolve().parents[2]
    try:
        recover(project, load_scheduler(project), apply=args.apply)
    except Exception as error:
        print("STOP: " + str(error), file=sys.stderr)
        print("If interrupted, rerun this same command; do not delete checkpoints or recovery receipts.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
