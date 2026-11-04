import argparse,json,os,subprocess,sys,time
from pathlib import Path
import storage
from common import META,unsafe,exclusive,recovery_ready
BUNDLE=Path(__file__).resolve().parent

def classify(rows,current=''):
    others=[r for r in rows if r[0]!=current]
    if any(r[1] in ('slide_next_gate','slide_next_setup') or r[1].startswith(('slide_inplace_','slide_expanded_','slide_roll_')) for r in others):
        raise RuntimeError('Another continuation is already active; no duplicate launch')
    return unsafe(others)

def submit(a,kind,dependency=None,delay=False):
    state=a.launch_state.resolve();state.mkdir(parents=True,exist_ok=True);(state/'logs').mkdir(exist_ok=True)
    pending=state/'pending_submission.json'
    if pending.exists():raise RuntimeError('Unresolved submission: '+str(pending))
    cmd=['sbatch','--parsable','--export=ALL','--exclude=midway2-0003,midway2-0011',
      '--output='+str(state/'logs'/(kind+'_%j.out')),'--error='+str(state/'logs'/(kind+'_%j.err'))]
    if dependency:cmd+=['--dependency=afterany:'+dependency]
    if delay:cmd+=['--begin=now+5minutes']
    common=[str(a.project.resolve()),str(a.recovery_state.resolve()),str(state)]
    cmd += [str(BUNDLE/('gate_job' if kind=='gate' else 'setup_job')),str(BUNDLE)]+(common if kind=='gate' else common[:2])
    storage.atomic_json(pending,dict(command=cmd,kind=kind,time=time.time()))
    job=subprocess.check_output(cmd,universal_newlines=True).strip().split(';')[0]
    if not job.isdigit():raise RuntimeError('Ambiguous sbatch response')
    storage.atomic_json(state/'last_submission.json',dict(job=job,kind=kind,time=time.time()))
    pending.unlink();storage.sync_dir(state);print('Submitted '+kind+' job '+job,flush=True)

def run(a):
    import fcntl
    storage.verify_bundle()
    state=a.launch_state.resolve();state.mkdir(parents=True,exist_ok=True)
    identity=dict(project=str(a.project.resolve()),recovery_state=str(a.recovery_state.resolve()),bundle=str(BUNDLE),target=40)
    with (state/'.launch.lock').open('a') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        p=state/'identity.json'
        if p.exists() and json.loads(p.read_text())!=identity:raise RuntimeError('Launch state belongs to another continuation')
        storage.atomic_json(p,identity)
        current=os.environ.get('SLURM_JOB_ID','') if a.mode=='gate' else ''
        rows=storage.queue();waiting=classify(rows,current)
        if a.mode=='submit':
            if (state/'last_submission.json').exists() and not a.resume:raise RuntimeError('Already submitted. Inspect logs; use --resume only after resolving a stopped run')
            if len(rows)>=300:raise RuntimeError('All 300 slots occupied. Retry when a slot becomes free; nothing submitted')
            submit(a,'gate',delay=True);return
        if len([r for r in rows if r[0]==current and r[1]=='slide_next_gate'])!=1:raise RuntimeError('Gate is not its active Slurm job')
        if waiting:
            print('Waiting for original recovery/005 jobs; leaving them unchanged: '+repr(waiting[:6]),flush=True)
            submit(a,'gate',dependency=current,delay=True);return
        project=a.project.resolve();recovery_ready(project,a.recovery_state.resolve())
        if not (project/META/'INITIALIZED.json').exists():
            with exclusive(project):
                storage.check_sources(project)
                subprocess.check_call([os.environ.get('RSCRIPT','Rscript'),'--vanilla',str(BUNDLE/'check_source.R'),str(BUNDLE),str(project)])
        submit(a,'setup',dependency=current)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['submit','gate'])
    for name in ('project','recovery_state','launch_state'):p.add_argument(name,type=Path)
    p.add_argument('--resume',action='store_true');a=p.parse_args()
    try:run(a)
    except Exception as e:print('STOP: '+str(e)+'. No continuation was started from an unvalidated checkpoint.',file=sys.stderr);sys.exit(1)
