"""One Slurm refill controller for slider EM and remaining simulation conditions.

Python >= 3.6, standard library only. No cancellation, package changes or fit math.
All submission state is durable; ambiguous submissions stop for inspection.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from contextlib import contextmanager

TERMINAL = set("COMPLETED CANCELLED FAILED TIMEOUT NODE_FAIL OUT_OF_MEMORY PREEMPTED BOOT_FAIL DEADLINE REVOKED".split())
EM_WORKERS = {"slide_parallel_fit", "slide_parallel_validate"}
SIM_WORKER = "slide_shared_sim"
SIM_SOURCES = ["script/sim/run_job.R", "script/sim/sim_workhorse.R", "script/sim/simulation_design.R",
               "script/sim/additive_init_genotypes.R", "script/scan_tissue_attempt/workhorse_utils.R",
               "script/scan_tissue_attempt/get_gene_annotations.R"]


def atomic(path, value, raw=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp." + str(os.getpid()))
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(value if raw else json.dumps(value, indent=2, sort_keys=True) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(str(temporary), str(path))


def read(path):
    with Path(path).open(encoding="utf-8") as stream:
        return json.load(stream)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(args):
    return subprocess.check_output([str(x) for x in args], universal_newlines=True).strip()


@contextmanager
def lock(path):
    import fcntl
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("a") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError("Another controller/migration holds {}. Retry in five minutes.".format(path))
        yield


def queue():
    rows = []
    raw = run(["squeue", "--noheader", "--array", "--states=all", "--user", os.environ["USER"],
               "--format=%i|%F|%200j|%T|%r"])
    for line in raw.splitlines():
        parts = [p.strip() for p in line.split("|")]
        if len(parts) != 5 or not re.fullmatch(r"[0-9]+(?:_[0-9]+)?", parts[0]) or not parts[1].isdigit():
            raise RuntimeError("Unrecognized expanded squeue record: " + line)
        if parts[3] not in TERMINAL:
            rows.append(dict(zip(("id", "array", "name", "state", "reason"), parts)))
    if len({r["id"] for r in rows}) != len(rows):
        raise RuntimeError("Duplicate queue records")
    return rows


def paths(project):
    return project / "results_slide_prior_em/shared_scheduler"


def plan_rows(project):
    plan = project / "script/sim/jobs_slide/submission_batches.csv"
    with plan.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    last = 0
    for i, row in enumerate(rows, 1):
        expected = "job/run_simulation_slide_batch_{}".format(i)
        if (int(row["batch"]) != i or int(row["array_start"]) != 0 or
                not 0 <= int(row["array_end"]) < 400 or int(row["offset"]) != last or
                int(row["first_job"]) != last + 1 or
                int(row["last_job"]) != last + int(row["array_end"]) + 1 or
                row["launcher"] != expected or not (project / expected).is_file()):
            raise RuntimeError("Invalid simulation submission plan")
        last = int(row["last_job"])
    if not rows:
        raise RuntimeError("Empty simulation plan")
    return rows


def bridge(batch):
    # Only future sbatch calls see this. Already submitted arrays keep Slurm's
    # saved script. A saved OLD launcher may still submit this harmless receipt
    # job, then submit the new launch_simulation_slide (which simply exits).
    return """#!/bin/bash
# Shared scheduler migration receipt, batch {batch}; no simulation array.
#SBATCH --job-name=susie_slide_bridge
#SBATCH --time=00:05:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256M
set -euo pipefail
project_dir="${{SUSIE_MIX_PROJECT_DIR:-/project2/mstephens/wdenault/susie_mix}}"
test -f "$project_dir/results_slide_prior_em/shared_scheduler/enabled.json"
echo 'Batch {batch} is owned by the shared EM/simulation controller. No duplicate array submitted.'
""".format(batch=batch)


def check_config(project, config):
    if config["schema"] != 1 or str(project) != config["project"]:
        raise RuntimeError("Shared scheduler configuration belongs to another project/version")
    for name, expected in config["source_hashes"].items():
        if digest(project / name) != expected:
            raise RuntimeError("Frozen simulation source changed: " + name)
    for b in config["bridges"]:
        name = "job/run_simulation_slide_batch_{}".format(b)
        if (project / name).read_text(encoding="utf-8") != bridge(b):
            raise RuntimeError("Migration bridge changed/regenerated: " + name)


def enable(project):
    root = paths(project)
    root.mkdir(parents=True, exist_ok=True)
    legacy = project / "simulation results/n250_r300_v1/launcher"
    migration_file = root / "migration.json"
    with lock(root / ".migration.lock"):
        if (root / "enabled.json").exists():
            config = read(root / "enabled.json")
            check_config(project, config)
            migration = read(migration_file)
            release_held(migration, migration_file)
            print("Sharing already enabled. Reuse the existing EM controller chain; do not submit another.")
            return
        rows = queue()
        controls = [r for r in rows if r["name"] == "slide_parallel_control"]
        if len(controls) != 1:
            raise RuntimeError("Expected exactly one existing EM controller (including its pending continuation); found {}. Inspect squeue first.".format(len(controls)))
        if controls[0]["state"] != "PENDING":
            raise RuntimeError("The EM controller is running its short check; retry enable in five minutes.")
        if any(r["name"].startswith("slide_recover_") or r["name"] in
               ("em_slide_prior_prepare", "em_susie_slide_prior") for r in rows):
            raise RuntimeError("Legacy EM chain still active")
        if migration_file.exists():
            migration = read(migration_file)
        else:
            controllers = [r for r in rows if r["name"] == "susie_slide_launch"]
            recorded = (legacy / "last_continuation_job_id.txt").read_text().strip() if (legacy / "last_continuation_job_id.txt").exists() else ""
            if len(controllers) > 1 or any(r["id"] != recorded or r["state"] != "PENDING" for r in controllers):
                raise RuntimeError("Simulation launcher is running or unrecognized; retry after its short submission step.")
            if any(r["reason"].startswith("JobHeld") for r in controllers + controls):
                raise RuntimeError("A launcher was already held; leave that user/admin hold unchanged.")
            migration = {"hold_ids": [r["id"] for r in controllers + controls], "released": [], "config": None}
            atomic(migration_file, migration)
        # Save ownership before hold so interrupted activation is recoverable.
        try:
            live_ids = {r["id"] for r in queue()}
            for job in migration["hold_ids"]:
                if job not in migration["released"] and job in live_ids:
                    run(["scontrol", "hold", job])
            live = queue()
            if any(r["name"] in ("susie_slide_launch", "slide_parallel_control") and
                   (r["state"] != "PENDING" or r["id"] not in migration["hold_ids"]) for r in live):
                raise RuntimeError("A launcher started/changed during handover; retry after its short check.")
            with lock(legacy / ".submit.lock"), lock(project / "results_slide_prior_em/.prepare.lock"):
                if migration["config"] is None:
                    batches = plan_rows(project)
                    last_batch = int((legacy / "last_batch.txt").read_text().strip())
                    last_array = (legacy / "last_array_job_id.txt").read_text().strip()
                    submissions = [line.split("\t") for line in (legacy / "submissions.tsv").read_text().splitlines()]
                    array_records = [(int(b), j) for b, kind, j in submissions if kind == "array"]
                    arrays = dict(array_records)
                    if (not 1 <= last_batch <= len(batches) or set(arrays) != set(range(1, last_batch + 1)) or
                            len(array_records) != len(arrays) or
                            arrays[last_batch] != last_array or any(not x.isdigit() for x in arrays.values())):
                        raise RuntimeError("Inconsistent legacy simulation submission records; no jobs changed")
                    if any(r["name"] == "susie_slide_sim" and r["array"] not in arrays.values() for r in live):
                        raise RuntimeError("Unregistered simulation array; inspect before adopting")
                    first = int(batches[last_batch - 1]["last_job"]) + 1
                    jobs = list(range(first, int(batches[-1]["last_job"]) + 1))
                    sources = SIM_SOURCES + ["script/sim/jobs_slide/sim_job_{}.R".format(j) for j in jobs]
                    datadir = os.environ.get("SUSIE_MIX_GTEX_DIR", "/project2/mstephens/gtex")
                    simulation_env = {"SUSIE_MIX_SIM_INPUT": os.environ.get("SUSIE_MIX_SIM_INPUT", "gtex"),
                                      "SUSIE_MIX_SIM_GENE": os.environ.get("SUSIE_MIX_SIM_GENE", ""),
                                      "SUSIE_MIX_GENOTYPE_DIR": os.environ.get("SUSIE_MIX_GENOTYPE_DIR", str(project / "temp_plink")),
                                      "SUSIE_MIX_GTEX_DIR": datadir,
                                      "SUSIE_MIX_PLINK": os.environ.get("SUSIE_MIX_PLINK", datadir + "/plink2"),
                                      "SUSIE_MIX_PLINK_MEMORY_MB": os.environ.get("SUSIE_MIX_PLINK_MEMORY_MB", "8000")}
                    if simulation_env["SUSIE_MIX_SIM_INPUT"] != "gtex" or simulation_env["SUSIE_MIX_SIM_GENE"]:
                        raise RuntimeError("This handover expects the submitted random-GTEx-gene simulation design (gtex mode, no fixed gene).")
                    config = {"schema": 1, "project": str(project), "jobs": jobs,
                              "legacy_arrays": list(arrays.values()), "legacy_last_batch": last_batch,
                              "simulation_env": simulation_env,
                              "bridges": list(range(last_batch + 1, len(batches) + 1)),
                              "source_hashes": {name: digest(project / name) for name in sources}}
                    # Snapshot originals before modifying any launcher on disk.
                    for b in config["bridges"]:
                        name = "run_simulation_slide_batch_{}".format(b)
                        atomic(root / "original_launchers" / name, (project / "job" / name).read_text(), raw=True)
                    migration["config"] = config
                    atomic(migration_file, migration)
                config = migration["config"]
                for b in config["bridges"]:
                    atomic(project / "job/run_simulation_slide_batch_{}".format(b), bridge(b), raw=True)
                check_config(project, config)
                atomic(root / "enabled.json", config)
        except Exception:
            # If any disk launchers were changed, keep the tiny pending launcher
            # held until enable is rerun to finish the recorded transaction.
            if migration["config"] is None:
                release_held(migration, migration_file)
                migration_file.unlink()
            else:
                print("Handover incomplete: rerun enable. Only recorded continuation launcher(s) remain held: " +
                      ", ".join(migration["hold_ids"]), file=sys.stderr)
            raise
        release_held(migration, migration_file)
        print("Sharing enabled: {} future simulation conditions registered; submitted arrays unchanged.".format(len(config["jobs"])))
        print("The existing EM controller chain adopts sharing on its next newly submitted continuation (normally one or two checks).")


def release_held(migration, path):
    live_ids = {r["id"] for r in queue()}
    for job in migration["hold_ids"]:
        if job not in migration["released"]:
            if job in live_ids:
                run(["scontrol", "release", job])
            migration["released"].append(job)
            atomic(path, migration)


def allocate(free, em_live, sim_live, sim_pending, limit, em_pending=True):
    """Fill toward an equal share, borrowing unused space; no existing preemption."""
    em = sim = 0
    for _ in range(max(0, free)):
        can_em = em_pending and em_live + em < limit
        can_sim = sim < sim_pending
        if can_em and (not can_sim or em_live + em <= sim_live + sim):
            em += 1
        elif can_sim:
            sim += 1
        else:
            break
    return em, sim


def sim_pending(project, config, state, rows):
    root = paths(project)
    active = {r["id"] for r in rows}
    owned = {a["array"] for a in state["arrays"]}
    if any(r["name"] == SIM_WORKER and r["array"] not in owned for r in rows):
        raise RuntimeError("Unregistered shared simulation worker; inspect pending submission")
    jobs = {j: [] for j in config["jobs"]}
    for a in state["arrays"]:
        for i, job in enumerate(a["jobs"], 1):
            jobs[job].append((a["array"] + "_" + str(i), a["batch"], i))
    pending = []
    for job, attempts in sorted(jobs.items()):
        owners = sum(task in active for task, _, _ in attempts)
        if owners > 1:
            raise RuntimeError("Overlapping live simulation assignments: " + str(job))
        if owners:
            continue
        done = root / "done/{}.json".format(job)
        if done.exists():
            receipt = read(done)
            if receipt["job"] != job or receipt["source_hash"] != config["source_hashes"]["script/sim/jobs_slide/sim_job_{}.R".format(job)]:
                raise RuntimeError("Invalid simulation completion receipt")
            continue
        if len(attempts) >= 3:
            raise RuntimeError("Simulation condition {} exhausted three attempts. Inspect shared_scheduler/batches; checkpoints retained.".format(job))
        pending.append(job)
    return pending


def submit(root, description, args):
    pending = root.parent / "parallel_pending_submission.txt"
    if pending.exists():
        raise RuntimeError("Unresolved submission: " + str(pending))
    atomic(pending, description, raw=True)
    try:
        answer = run(["sbatch", "--parsable"] + args)
    except subprocess.CalledProcessError:
        pending.unlink()
        raise
    job = answer.split(";")[0]
    if not job.isdigit():
        raise RuntimeError("Ambiguous sbatch response; inspect " + str(pending))
    return job  # Caller persists registration BEFORE clearing the marker.


def clear_pending(root):
    (root.parent / "parallel_pending_submission.txt").unlink()


def control(project):
    root = paths(project)
    with lock(root.parent / ".prepare.lock"):
        if (root.parent / "parallel_pending_submission.txt").exists():
            raise RuntimeError("Unresolved submission; inspect parallel_pending_submission.txt")
        config = read(root / "enabled.json")
        check_config(project, config)
        current = os.environ["SLURM_JOB_ID"]
        rows = queue()
        if len([r for r in rows if r["id"] == current and r["name"] == "slide_parallel_control"]) != 1:
            raise RuntimeError("Queue does not include this controller")
        if any((r["name"] == "slide_parallel_control" and r["id"] != current) or
               r["name"].startswith("slide_recover_") or r["name"] in
               ("em_slide_prior_prepare", "em_susie_slide_prior") for r in rows):
            raise RuntimeError("Another EM controller or legacy worker is active")
        if any(r["name"] == "susie_slide_sim" and r["array"] not in config["legacy_arrays"] for r in rows):
            raise RuntimeError("An unregistered legacy simulation array appeared after handover; inspect before refilling")
        state_file = root / "state.json"
        state = read(state_file) if state_file.exists() else {"arrays": []}
        pending = sim_pending(project, config, state, rows)
        em_live = sum(r["name"] in EM_WORKERS for r in rows)
        sim_live = sum(r["name"] in (SIM_WORKER, "susie_slide_sim") for r in rows)
        legacy = any(r["name"] in ("susie_slide_launch", "susie_slide_bridge") for r in rows)
        # While the OLD saved simulation launcher still exists, reserve its
        # bridge + next-launcher submissions as well. They will quickly drain.
        free = max(0, 300 - len(rows) - 1 - (2 if legacy else 0))
        limit = int(os.environ.get("SUSIE_SLIDE_WORKERS", "298"))
        poll = int(os.environ.get("SUSIE_SLIDE_POLL_MINUTES", "5"))
        if not 1 <= limit <= 298 or not 1 <= poll <= 60:
            raise RuntimeError("Invalid worker limit or polling interval")
        em_capacity, _ = allocate(free, em_live, sim_live, len(pending), limit)
        qfile = root / "queue.txt"
        atomic(qfile, "".join("{id}|{array}|{name}\n".format(**r) for r in rows if r["id"] != current), raw=True)
        cli = project / "script/scan_tissue_attempt/run_slide_prior_parallel.R"
        response = run(["Rscript", "--vanilla", cli, "plan", project, "latest", str(em_capacity), qfile]).splitlines()
        if len(response) < 2 or response[0] not in ("batch", "waiting", "complete") or not response[1].isdigit():
            raise RuntimeError("Invalid EM plan response")
        extra = []
        exclude = os.environ.get("SUSIE_SLIDE_EXCLUDE_NODES", "")
        if exclude:
            if not re.fullmatch(r"[A-Za-z0-9_.-]+(?:,[A-Za-z0-9_.-]+)*", exclude):
                raise RuntimeError("Invalid node exclusions")
            extra = ["--exclude=" + exclude]
        common = extra + ["--export=ALL", "--chdir=" + str(project / "job")]
        used_em = 0
        if response[0] == "batch":
            if len(response) != 5 or response[4] not in ("fit", "validate") or not response[3].isdigit():
                raise RuntimeError("Invalid EM worker plan")
            batch, tasks, action = Path(response[2]), int(response[3]), response[4]
            if not batch.is_dir() or not 0 < tasks <= em_capacity:
                raise RuntimeError("EM worker plan exceeds capacity")
            job = submit(root, str(batch), common + ["--array=1-" + str(tasks), "--job-name=slide_parallel_" + action,
                         "--output=" + str(batch / "worker_%A_%a.out"), "--error=" + str(batch / "worker_%A_%a.err"),
                         str(project / "job/slide_prior_parallel_worker"), str(project), response[1], str(batch)])
            atomic(batch / "array_job_id.txt", job + "\n", raw=True)
            clear_pending(root)
            used_em = tasks
            print("Submitted EM {} array {}: {} tasks".format(action, job, tasks))
        # Lend the simulations capacity that EM could not use (including when
        # the current iteration is waiting on its last genes, or target reached).
        sim_count = min(len(pending), free - used_em)
        if sim_count:
            batch = root / "batches" / (str(len(state["arrays"]) + 1).zfill(6))
            batch.mkdir(parents=True, exist_ok=True)
            assignment = {"jobs": pending[:sim_count]}
            atomic(batch / "assignment.json", assignment)
            job = submit(root, str(batch), common + ["--array=1-" + str(sim_count),
                         "--output=" + str(batch / "worker_%A_%a.out"), "--error=" + str(batch / "worker_%A_%a.err"),
                         str(project / "job/slide_shared_sim_worker"), str(project), str(batch)])
            state["arrays"].append({"array": job, "batch": str(batch), "jobs": assignment["jobs"]})
            atomic(state_file, state)
            clear_pending(root)
            print("Submitted simulation array {}: {} conditions (existing checkpoints reused)".format(job, sim_count))
        print("Queue snapshot: {} EM workers, {} simulation workers; {} additions; 300-job total budget.".format(em_live, sim_live, used_em + sim_count))
        if response[0] == "complete" and not pending and sim_live == 0 and not legacy:
            print("EM target and adopted simulations complete. No continuation submitted.")
            return
        job = submit(root, "continuation after " + current, common + ["--dependency=afterany:" + current,
                     "--begin=now+{}minutes".format(poll),
                     "--output=" + str(project / "job/slide_parallel_control_%j.out"),
                     "--error=" + str(project / "job/slide_parallel_control_%j.err"),
                     str(project / "job/resume_slide_prior_em"), "latest"])
        atomic(root.parent / "parallel_continuation_job_id.txt", job + "\n", raw=True)
        clear_pending(root)
        print("Queued shared controller {}; eligible in {} minutes.".format(job, poll))


def worker(project, batch, task):
    root = paths(project)
    config = read(root / "enabled.json")
    assignment = read(batch / "assignment.json")["jobs"]
    if task < 1 or task > len(assignment):
        raise RuntimeError("Invalid simulation task index")
    job = assignment[task - 1]
    if job not in config["jobs"]:
        raise RuntimeError("Simulation condition was not adopted")
    name = "script/sim/jobs_slide/sim_job_{}.R".format(job)
    for source in [name] + SIM_SOURCES:
        if digest(project / source) != config["source_hashes"][source]:
            raise RuntimeError("Frozen simulation source changed: " + source)
    with lock(root / "locks/{}.lock".format(job)):
        atomic(batch / "started_{}.json".format(task), {"job": job, "slurm": os.environ["SLURM_JOB_ID"]})
        # R itself validates checkpoints and skips completed replicates.
        environment = dict(os.environ)
        environment.update(config["simulation_env"])
        code = subprocess.call(["Rscript", "--vanilla", str(project / name)], env=environment)
        atomic(batch / "exit_{}.json".format(task), {"job": job, "exit_code": code})
        if code:
            raise RuntimeError("Simulation condition {} exited {}".format(job, code))
        atomic(root / "done/{}.json".format(job), {"job": job, "source_hash": config["source_hashes"][name]})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("enable", "control", "worker"))
    parser.add_argument("project", type=Path)
    parser.add_argument("batch", nargs="?", type=Path)
    parser.add_argument("task", nargs="?", type=int)
    args = parser.parse_args()
    project = args.project.resolve()
    try:
        if args.mode == "enable":
            enable(project)
        elif args.mode == "control":
            control(project)
        else:
            worker(project, args.batch.resolve(), args.task)
    except Exception as error:
        print("STOP: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
