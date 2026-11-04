"""Offline, resumable slider checkpoint maintenance. Python 3.6+, stdlib.

All mutations require no live slider jobs and an exclusive preparation lock.
Original bytes remain in verified scratch escrow until replacement is committed.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import threading
import time
from contextlib import contextmanager

BUNDLE = Path(__file__).resolve().parent
TERMINAL = set("COMPLETED CANCELLED FAILED TIMEOUT NODE_FAIL OUT_OF_MEMORY PREEMPTED BOOT_FAIL DEADLINE REVOKED".split())

def sync_dir(path):
    if os.name != "nt":
        fd = os.open(str(path), os.O_RDONLY)
        try: os.fsync(fd)
        finally: os.close(fd)

def atomic_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", encoding="utf-8") as f:
        json.dump(value, f, indent=2, sort_keys=True); f.write("\n"); f.flush(); os.fsync(f.fileno())
    os.replace(str(temp), str(path)); sync_dir(path.parent)

def hashes(path):
    sha, md = hashlib.sha256(), hashlib.md5()
    with Path(path).open("rb") as f:
        while True:
            b = f.read(8*1024*1024)
            if not b: break
            sha.update(b); md.update(b)
    return {"sha256": sha.hexdigest(), "md5": md.hexdigest()}

def stamp(path):
    st = Path(path).stat()
    return [st.st_size, st.st_mtime_ns]

def copy_verified(source, destination, expected=None):
    destination = Path(destination)
    with Path(source).open("rb") as src, destination.open("wb") as dst:
        shutil.copyfileobj(src, dst, 8*1024*1024); dst.flush(); os.fsync(dst.fileno())
    expected = expected or hashes(source)
    if hashes(destination) != expected: raise RuntimeError("Copy verification failed: " + str(destination))
    sync_dir(destination.parent)
    return expected

def queue():
    raw = subprocess.check_output(["squeue", "--noheader", "--array", "--states=all", "--user", os.environ["USER"],
                                   "--format=%i|%200j|%T"], universal_newlines=True)
    rows = []
    for line in raw.splitlines():
        fields = [x.strip() for x in line.split("|")]
        if len(fields) != 3: raise RuntimeError("Unrecognized squeue output")
        if fields[2] not in TERMINAL: rows.append(fields)
    return rows

@contextmanager
def maintenance_lock(project):
    import fcntl
    lock = project / "results_slide_prior_em/.prepare.lock"
    with lock.open("a") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        unsafe = [r for r in queue() if r[1].startswith(("slide_parallel_", "slide_recover_", "slide_expanded_"))
                  or r[1] in ("em_slide_prior_prepare", "em_susie_slide_prior", "slide_prior_em")]
        if unsafe: raise RuntimeError("Slider jobs still active; let them finish. No jobs cancelled: " + repr(unsafe[:10]))
        yield

def run_r(*args):
    subprocess.check_call([os.environ.get("RSCRIPT", "Rscript"), "--vanilla", str(BUNDLE / "maintenance.R")] + [str(x) for x in args])

def verify_bundle():
    for line in (BUNDLE / "SHA256SUMS").read_text().splitlines():
        digest, relative = line.split("  ", 1)
        if hashes(BUNDLE / relative)["sha256"] != digest: raise RuntimeError("Bundle checksum mismatch: " + relative)

def check_sources(project):
    expected = json.loads((BUNDLE / "expected_sources.json").read_text())
    for rel, digest in expected.items():
        text = (project / rel).read_text().replace("\r\n", "\n")
        current = hashlib.sha256(text.encode()).hexdigest()
        patched = hashes(BUNDLE / "patch" / rel)["sha256"]
        if current not in (digest, patched): raise RuntimeError("Unrecognized deployed source; retained unchanged: " + rel)

def install(project, state):
    check_sources(project)
    rel = "script/scan_tissue_attempt/workhorse_slide_prior_em.R"
    source, dest = BUNDLE / "patch" / rel, project / rel
    backup = state / "original_workhorse.R"
    if not backup.exists(): copy_verified(dest, backup)
    if hashes(source) != hashes(dest):
        temp = dest.with_name(dest.name + ".storage-tmp")
        copy_verified(source, temp)
        os.replace(str(temp), str(dest)); sync_dir(dest.parent)

class Ledger:
    def __init__(self, path):
        self.path = path
        self.mutex = threading.Lock()
        self.db = sqlite3.connect(str(path), check_same_thread=False)
        self.db.execute("PRAGMA synchronous=FULL")
        self.db.execute("CREATE TABLE IF NOT EXISTS files (relative TEXT PRIMARY KEY, receipt TEXT NOT NULL)")
        self.db.commit()
    def get(self, relative):
        with self.mutex:
            row = self.db.execute("SELECT receipt FROM files WHERE relative=?", (relative,)).fetchone()
        return json.loads(row[0]) if row else None
    def put(self, relative, receipt):
        with self.mutex:
            self.db.execute("INSERT OR REPLACE INTO files VALUES (?,?)", (relative, json.dumps(receipt)))
            self.db.commit()
    def export(self, path):
        keys = ["relative", "old_md5", "new_md5", "old_bytes", "new_bytes", "removed_calls"]
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
            for rel, value in self.db.execute("SELECT relative, receipt FROM files ORDER BY relative"):
                r = json.loads(value); w.writerow(dict(relative=rel, **{k:r[k] for k in keys[1:]}))
            f.flush(); os.fsync(f.fileno())

def contained(path, root):
    try: path.resolve().relative_to(root.resolve()); return True
    except ValueError: return False

def transact(project, relative, stage, ledger, compact=None):
    """One verified replacement. The same stage is resumed after an interruption."""
    src = project / relative
    if not contained(src, project) or src.is_symlink(): raise RuntimeError("Unsafe/symlink result path: " + relative)
    done = ledger.get(relative)
    stage.mkdir(parents=True, exist_ok=True)
    journal = stage / "transaction.json"
    if done:
        if not src.exists() or stamp(src) != done["stamp"]: raise RuntimeError("Committed result changed: " + relative)
        # A crash after commit can leave only escrow cleanup outstanding.
        for name in ("original.rds", "compact.rds", "report.csv", "transaction.json"):
            p = stage / name
            if p.exists(): p.unlink()
        return done
    if journal.exists():
        t = json.loads(journal.read_text())
        if t["relative"] != relative: raise RuntimeError("Scratch transaction belongs to another file")
    else:
        if not src.is_file(): raise RuntimeError("Result disappeared without an escrow transaction: " + relative)
        t = {"relative":relative, "original_stamp":stamp(src), "old_hash":hashes(src), "old_bytes":src.stat().st_size,
             "stage":"start"}
        atomic_json(journal,t)
    original, packed = stage / "original.rds", stage / "compact.rds"
    if t["stage"] == "start":
        if stamp(src) != t["original_stamp"] or hashes(src) != t["old_hash"]: raise RuntimeError("Original changed during maintenance")
        copy_verified(src, original, t["old_hash"])
        if compact is None:
            log = stage / "R.log"
            with log.open("w") as out:
                code = subprocess.call([os.environ.get("RSCRIPT","Rscript"),"--vanilla",str(BUNDLE/"maintenance.R"),
                                        "file",str(original),str(packed),str(stage/"report.csv")],stdout=out,stderr=subprocess.STDOUT)
            if code: raise RuntimeError("Compaction failed; original retained. See " + str(log))
        else: compact(original, packed, stage/"report.csv")
        with (stage/"report.csv").open() as f: t["removed_calls"] = int(next(csv.DictReader(f))["removed_calls"])
        # A no-op needs no rewrite; gzip may differ across R versions.
        if not t["removed_calls"]:
            t["new_hash"] = t["old_hash"]; t["new_bytes"] = t["old_bytes"]
            t["stage"] = "unchanged"
        else:
            with packed.open("r+b") as f: os.fsync(f.fileno())
            t["new_hash"] = hashes(packed); t["new_bytes"] = packed.stat().st_size
            if t["new_bytes"] >= t["old_bytes"]: raise RuntimeError("Compaction did not reduce disk size; original retained")
            t["stage"] = "prepared"
        atomic_json(journal,t)
    if t["stage"] in ("prepared", "vacated"):
        # Never discard the original until BOTH escrow files have verified hashes.
        if hashes(original) != t["old_hash"] or hashes(packed) != t["new_hash"]:
            raise RuntimeError("Escrow verification failed; do not remove scratch staging")
        if src.exists():
            current = hashes(src)
            if current == t["new_hash"]:
                t["stage"] = "replaced"
            elif current == t["old_hash"] and stamp(src) == t["original_stamp"]:
                if src.stat().st_nlink != 1: raise RuntimeError("Hard-linked source cannot release quota safely")
                src.unlink(); sync_dir(src.parent)
                t["stage"] = "vacated"; atomic_json(journal,t)
            else: raise RuntimeError("Unexpected result at transaction target")
        if t["stage"] != "replaced":
            temp = src.with_name(".storage_compact_" + src.name)
            copy_verified(packed,temp,t["new_hash"])
            os.replace(str(temp),str(src)); sync_dir(src.parent)
            t["stage"] = "replaced"; atomic_json(journal,t)
    if t["stage"] not in ("replaced", "unchanged"): raise RuntimeError("Unknown transaction stage")
    if hashes(src) != t["new_hash"]: raise RuntimeError("Final result verification failed")
    receipt = dict(old_md5=t["old_hash"]["md5"],new_md5=t["new_hash"]["md5"],old_bytes=t["old_bytes"],
                   new_bytes=t["new_bytes"],removed_calls=t["removed_calls"],stamp=stamp(src))
    ledger.put(relative,receipt)  # Durable receipt precedes removal of the original escrow.
    for p in (original,packed,stage/"report.csv",journal):
        if p.exists(): p.unlink()
    return receipt

def inventory(project, through=None):
    files=[]
    root=project/"results_slide_prior_em"
    for d in sorted(root.glob("iteration_[0-9][0-9][0-9]*")):
        iteration=int(d.name.split("_")[1])
        if through is not None and iteration > through: continue
        if (d/"storage_retired.json").exists(): continue
        for p in sorted((d/"results").glob("*.rds")):
            if p.is_symlink() or not contained(p,project): raise RuntimeError("Result is a symlink or escapes project")
            files.append(p.relative_to(project).as_posix())
    return files

def setup_state(project,state,kind,extra):
    state.mkdir(parents=True,exist_ok=True)
    expected=dict(project=str(project),kind=kind,**extra)
    p=state/"identity.json"
    if p.exists():
        if json.loads(p.read_text()) != expected: raise RuntimeError("State belongs to another operation")
    else: atomic_json(p,expected)

def compact_all(args):
    project,state,scratch=args.project.resolve(),args.state.resolve(),args.scratch.resolve()
    if contained(state,project) or contained(scratch,project): raise RuntimeError("State and scratch must be outside the project")
    verify_bundle(); check_sources(project)
    setup_state(project,state,"compact",dict(expected=args.expected,scratch=str(scratch)))
    if (state/"COMPLETE.json").exists():
        print("Compaction already complete. See",state/"validation_report.csv"); return
    with maintenance_lock(project):
        if not (state/"snapshot.rds").exists(): run_r("snapshot",BUNDLE,project,state,args.expected)
        planfile=state/"files.json"
        if planfile.exists(): files=json.loads(planfile.read_text())
        else: files=inventory(project); atomic_json(planfile,files)
        scratch.mkdir(parents=True,exist_ok=True)
        owner=dict(state=str(state),project=str(project))
        if (scratch/"owner.json").exists() and json.loads((scratch/"owner.json").read_text()) != owner:
            raise RuntimeError("Scratch belongs to another maintenance operation")
        atomic_json(scratch/"owner.json",owner)
        ledger=Ledger(state/"receipts.sqlite")
        # Per-file deterministic stages survive Slurm wall-time interruption.
        def one(rel):
            stage=scratch/hashlib.sha256(rel.encode()).hexdigest()[:24]
            try: return transact(project,rel,stage,ledger)
            except Exception as e: raise RuntimeError(rel+": "+str(e))
        n=0; before=after=0
        # Bound outstanding tasks so a failure does not start thousands more.
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            it=iter(files); pending={}
            for _ in range(args.workers):
                rel=next(it,None)
                if rel is not None: pending[pool.submit(one,rel)]=rel
            while pending:
                done,_=concurrent.futures.wait(pending,return_when=concurrent.futures.FIRST_COMPLETED)
                for future in done:
                    pending.pop(future); r=future.result(); n+=1; before+=r["old_bytes"]; after+=r["new_bytes"]
                    if n%25==0 or n==len(files): print("Verified {}/{}; {:.2f} GiB released".format(n,len(files),(before-after)/2**30),flush=True)
                    rel=next(it,None)
                    if rel is not None: pending[pool.submit(one,rel)]=rel
        # Recheck exactly the original result inventory: no silent external additions.
        if set(inventory(project)) != set(files): raise RuntimeError("Result inventory changed")
        ledger.export(state/"receipts.csv")
        install(project,state)
        run_r("finalize",BUNDLE,project,state)
        atomic_json(state/"COMPLETE.json",dict(project=str(project),files=n,old_bytes=before,new_bytes=after,
            saved_bytes=before-after,latest=args.expected,finished=time.time()))
        print("COMPACTION COMPLETE. Scientific fields round-trip identical; priors unchanged.",flush=True)

def retire(args):
    project,state=args.project.resolve(),args.state.resolve()
    if contained(state,project): raise RuntimeError("Retirement state must be outside project")
    files=inventory(project,args.through)
    total=sum((project/p).stat().st_size for p in files)
    print("Retire {} full-result files through iteration {:03d}: {:.2f} TiB. Keep all priors/counts/summaries.".format(len(files),args.through,total/2**40))
    if not args.apply:
        print("PLAN ONLY. Add --apply to execute this irreversible removal of historical full posteriors."); return
    verify_bundle(); check_sources(project)
    setup_state(project,state,"retire",dict(through=args.through,expected=args.expected))
    if (state/"COMPLETE.json").exists(): print("Retirement already complete"); return
    with maintenance_lock(project):
        if not (state/"retire_plan.rds").exists(): run_r("retire-prepare",BUNDLE,project,state,args.through,args.expected)
        pf=state/"files.json"
        if pf.exists(): files=json.loads(pf.read_text())
        else: atomic_json(pf,files)
        dependencies=project/"results_slide_prior_em/storage_dependencies.json"
        if dependencies.exists():
            deps=json.loads(dependencies.read_text())
            if any(int(x["source_iteration"])<=args.through for x in deps): raise RuntimeError("An expanded branch still depends on a selected source checkpoint")
        for rel in files:
            p=project/rel
            if not contained(p,project) or p.is_symlink() or p.parent.name!="results": raise RuntimeError("Unsafe retirement target")
            if p.exists(): p.unlink()
        run_r("retire-finish",BUNDLE,project,state)
        for i in range(args.through+1):
            d=project/"results_slide_prior_em"/"iteration_{:03d}".format(i)
            cert=d/"recovery/parallel/COMPLETE.rds"
            if cert.exists(): os.replace(str(cert),str(cert.with_name("COMPLETE.before_storage_retirement.rds")))
            atomic_json(d/"storage_retired.json",dict(state=str(state),through=args.through,retired=time.time()))
        atomic_json(state/"COMPLETE.json",dict(files=len(files),through=args.through,project=str(project)))
        print("Retirement complete. Prior history and retained checkpoint validation preserved.")

def resume(args):
    project,state=args.project.resolve(),args.state.resolve()
    completed=json.loads((state/"COMPLETE.json").read_text())
    if completed.get("project")!=str(project) or "saved_bytes" not in completed: raise RuntimeError("Compaction has not completed")
    with maintenance_lock(project):
        check_sources(project)
        env=os.environ.copy()
        env.update(SUSIE_MIX_PROJECT_DIR=str(project),SUSIE_SLIDE_TARGET_ITERATION=str(args.target),
                   SUSIE_SLIDE_WORKERS="298",SUSIE_SLIDE_POLL_MINUTES="5",SUSIE_SLIDE_MAX_ITER="5000",
                   SUSIE_SLIDE_RECORD_ELBO="0",SUSIE_SLIDE_EXCLUDE_NODES="midway2-0003,midway2-0011")
        dependency=["--dependency=afterany:"+os.environ["SLURM_JOB_ID"]] if os.environ.get("SLURM_JOB_ID") else ["--begin=now+1minutes"]
        subprocess.check_call(["sbatch"]+dependency+["--exclude=midway2-0003,midway2-0011","--chdir="+str(project/"job"),
            "--output="+str(project/"job/slide_parallel_control_%j.out"),
            "--error="+str(project/"job/slide_parallel_control_%j.err"),str(project/"job/resume_slide_prior_em"),"latest"],env=env)

def main():
    p=argparse.ArgumentParser(description=__doc__); subs=p.add_subparsers(dest="mode")
    for name in ("compact","retire","resume"):
        q=subs.add_parser(name); q.add_argument("project",type=Path); q.add_argument("state",type=Path)
        if name=="compact":
            q.add_argument("scratch",type=Path);q.add_argument("--workers",type=int,default=4);q.add_argument("--expected",type=int,default=5)
        elif name=="retire":
            q.add_argument("--through",type=int,required=True);q.add_argument("--expected",type=int,default=5);q.add_argument("--apply",action="store_true")
        else: q.add_argument("--target",type=int,default=5)
    a=p.parse_args()
    if a.mode=="compact":
        if not 1<=a.workers<=16: p.error("workers must be 1..16")
        compact_all(a)
    elif a.mode=="retire": retire(a)
    elif a.mode=="resume": resume(a)
    else: p.error("Choose compact, retire or resume")
if __name__=="__main__":
    try: main()
    except Exception as e:
        print("STOP: "+str(e)+"\nKeep the state and scratch directories. Rerun the same command to resume.",file=sys.stderr)
        sys.exit(1)
