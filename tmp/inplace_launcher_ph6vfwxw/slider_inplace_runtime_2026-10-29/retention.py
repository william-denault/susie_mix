"""Controller-only, restartable retirement of old full slider posteriors."""
import argparse, csv, json, os, stat, subprocess, sys, time
from pathlib import Path
from contextlib import contextmanager
import storage
from common import RUNTIME,CLEANUP,unsafe,check_dependencies

def rrun(project, mode, state, current=None):
    args=[os.environ.get('RSCRIPT','Rscript'),'--vanilla',str(project/RUNTIME/'retention.R'),mode,str(project),str(state)]
    if current is not None: args.append(str(current))
    subprocess.check_call(args)

def read_csv(path):
    with path.open(newline='',encoding='utf-8') as f: return list(csv.DictReader(f))

def allowed_dirs(project,current):
    dirs=[project/'results_slide_prior_em'/('iteration_%03d'%i) for i in range(current-1)]
    return dirs,project

def safe_directory(d,allowed):
    if d not in allowed or d.is_symlink() or d.resolve()!=d or not d.is_dir():
        raise RuntimeError('Unapproved retirement directory: '+str(d))
    result=d/'results'
    if result.is_symlink() or result.resolve()!=result or not result.is_dir():
        raise RuntimeError('Results directory is missing or redirected: '+str(result))

def plan_files(project,state,current):
    allowed,source=allowed_dirs(project,current)
    dirs=[Path(r['directory']) for r in read_csv(state/'directories.csv')]
    if dirs!=allowed: raise RuntimeError('Retirement scope differs from the two-checkpoint policy')
    entries=[]
    for d in dirs:
        safe_directory(d,allowed)
        for p in sorted((d/'results').iterdir()):
            s=p.lstat()
            if p.suffix!='.rds' or not stat.S_ISREG(s.st_mode) or p.resolve()!=p:
                raise RuntimeError('Unexpected object in results; retained: '+str(p))
            entries.append(dict(path=str(p),size=s.st_size,mtime_ns=s.st_mtime_ns,device=s.st_dev,inode=s.st_ino))
    inputs={p.name:storage.hashes(p)['sha256'] for p in (state/'snapshot.rds',state/'directories.csv')}
    plan=dict(project=str(project),current=current,directories=[str(d) for d in dirs],files=entries,inputs=inputs)
    storage.atomic_json(state/'plan.json',plan)
    return plan

def verify_plan(project,state,current):
    plan=json.loads((state/'plan.json').read_text())
    allowed,_=allowed_dirs(project,current)
    if plan['project']!=str(project) or plan['current']!=current or [Path(d) for d in plan['directories']]!=allowed:
        raise RuntimeError('Stored retirement plan belongs to another run')
    for name,digest in plan['inputs'].items():
        if name not in ('snapshot.rds','directories.csv') or storage.hashes(state/name)['sha256']!=digest:
            raise RuntimeError('Retirement metadata changed')
    for d in allowed: safe_directory(d,allowed)
    return plan

def execute_plan(project,state,current,after_unlink=None):
    plan=verify_plan(project,state,current)
    allowed={str(Path(d)/'results') for d in plan['directories']}
    expected={r['path'] for r in plan['files']}
    if len(expected)!=len(plan['files']): raise RuntimeError('Duplicate retirement entries')
    # Check the WHOLE remaining scope before the first deletion on every retry.
    for directory in allowed:
        for p in Path(directory).iterdir():
            if str(p) not in expected: raise RuntimeError('Unplanned file appeared; retained: '+str(p))
    for row in plan['files']:
        p=Path(row['path'])
        if str(p.parent) not in allowed or p.suffix!='.rds': raise RuntimeError('Invalid deletion target')
        if p.exists() or p.is_symlink():
            s=p.lstat()
            if not stat.S_ISREG(s.st_mode) or p.resolve()!=p or (s.st_size,s.st_mtime_ns,s.st_dev,s.st_ino)!=(row['size'],row['mtime_ns'],row['device'],row['inode']):
                raise RuntimeError('Historical file changed; retained: '+str(p))
    count=0
    for row in plan['files']:
        p=Path(row['path'])
        if p.exists():
            # Recheck immediately before unlink, even though both controllers
            # are excluded by locks and the live-queue guard.
            s=p.lstat()
            if not stat.S_ISREG(s.st_mode) or p.resolve()!=p or (s.st_size,s.st_mtime_ns,s.st_dev,s.st_ino)!=(row['size'],row['mtime_ns'],row['device'],row['inode']):
                raise RuntimeError('Historical file changed before unlink: '+str(p))
            p.unlink(); count+=1
            if after_unlink: after_unlink(p)
            if count%1000==0: print('Retired %d additional files before iteration %03d'%(count,current),flush=True)
    for directory in allowed: storage.sync_dir(directory)
    return plan

def archive_certificates(plan):
    for name in plan['directories']:
        d=Path(name); cert=d/'recovery/parallel/COMPLETE.rds'
        archived=cert.with_name('COMPLETE.before_rolling_retirement.rds')
        if cert.exists():
            if archived.exists():
                if storage.hashes(cert)!=storage.hashes(archived): raise RuntimeError('Conflicting archived certificate')
                cert.unlink()
            else: os.replace(str(cert),str(archived))
        storage.atomic_json(d/'storage_retired.json',dict(policy='two full checkpoints',before_iteration=plan['current'],
             history_preserved=True,full_results_retired=True))

def check_queue():
    current=os.environ.get('SLURM_JOB_ID','')
    rows=storage.queue()
    if not current or len([r for r in rows if r[0]==current and r[1]=='slide_inplace_control'])!=1:
        raise RuntimeError('Cleanup must run inside the active in-place controller')
    live=unsafe(rows,current)
    if live: raise RuntimeError('Slider jobs are still active; no results removed: '+repr(live[:8]))

@contextmanager
def barrier_lock(project,current):
    import fcntl
    # The shell controller MUST own the preparation lock throughout this call.
    with (project/'results_slide_prior_em/.prepare.lock').open('a') as probe:
        try: fcntl.flock(probe.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError: pass
        else: raise RuntimeError('Missing controller preparation lock')
    with (project/'results_slide_prior_em/.rolling.lock').open('a') as own:
        fcntl.flock(own.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        check_queue()
        check_dependencies(project)
        yield

def step(project,current,queue_file):
    project=project.resolve()
    if not 6<=current<=40: raise RuntimeError('Iteration must be 006..040')
    subprocess.check_call([sys.executable,str(project/RUNTIME/'verify_inplace.py'),str(project)])
    state=project/CLEANUP/('before_%03d'%current)
    if (state/'COMPLETE.json').exists(): return
    with barrier_lock(project,current):
        state.mkdir(parents=True,exist_ok=True)
        if not (state/'plan.json').exists():
            rrun(project,'prepare',state,current)
            plan_files(project,state,current)
        verify_plan(project,state,current)
        rrun(project,'verify',state)
        check_queue()
        plan=execute_plan(project,state,current)
        rrun(project,'finish',state)
        archive_certificates(plan)
        storage.atomic_json(state/'COMPLETE.json',dict(iteration=current,keep=[current-1,current],
             retired_files=len(plan['files']),unlinked_logical_bytes=sum(r['size'] for r in plan['files']),
             prior_history_preserved=True,finished=time.time()))
        print('Cleanup complete; retaining full iterations %03d and %03d.'%(current-1,current),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['step'])
    p.add_argument('project',type=Path);p.add_argument('current',type=int);p.add_argument('queue_file',type=Path)
    a=p.parse_args()
    try: step(a.project,a.current,a.queue_file)
    except Exception as e:
        print('STOP: '+str(e)+'; cleanup is restartable, and no new workers were submitted.',file=sys.stderr);sys.exit(1)
