"""Queue a read-only gate now; start expansion only after original 005 finishes."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
import storage

BUNDLE=Path(__file__).resolve().parent
EXCLUDE='midway2-0003,midway2-0011'

def is_slider(name):
    return name.startswith(('slide_parallel_','slide_expanded_','slide_recover_','slide_storage_')) or name in (
        'em_slide_prior_prepare','em_susie_slide_prior','slide_prior_em','slide_expand_setup','slide_roll_setup')

def classify(rows,current=''):
    others=[r for r in rows if r[0]!=current]
    if any(r[1]=='slide_roll_gate' or r[1].startswith('slide_expanded_') or r[1]=='slide_roll_setup' for r in others):
        raise RuntimeError('Another expanded run/gate is active; no duplicate submission')
    return [r for r in others if is_slider(r[1])]

def submit_job(a,kind,dependency=None,delay=False):
    state=a.launch_state.resolve(); state.mkdir(parents=True,exist_ok=True)
    (state/'logs').mkdir(exist_ok=True)
    pending=state/'pending_submission.json'
    if pending.exists(): raise RuntimeError('Unresolved submission; inspect '+str(pending)+' and squeue before retrying')
    cmd=['sbatch','--parsable','--exclude='+EXCLUDE,'--export=ALL',
         '--output='+str(state/'logs'/(kind+'_%j.out')),'--error='+str(state/'logs'/(kind+'_%j.err'))]
    if dependency: cmd+=['--dependency=afterany:'+dependency]
    if delay: cmd+=['--begin=now+5minutes']
    common=[str(a.project.resolve()),str(a.storage_state.resolve()),str(a.branch.resolve()),str(state)]
    if kind=='gate': cmd += [str(BUNDLE/'gate_job'),str(BUNDLE)]+common
    else: cmd += [str(BUNDLE/'expansion_job'),str(BUNDLE)]+common[:3]
    storage.atomic_json(pending,dict(command=cmd,kind=kind,time=time.time()))
    # Ambiguous/error responses leave the marker: never guess and resubmit.
    result=subprocess.check_output(cmd,universal_newlines=True).strip().split(';')[0]
    if not result.isdigit(): raise RuntimeError('Ambiguous sbatch response')
    storage.atomic_json(state/'last_submission.json',dict(job=result,kind=kind,time=time.time()))
    pending.unlink(); storage.sync_dir(state)
    print('Submitted '+kind+' job '+result,flush=True)

def run(a):
    import fcntl
    storage.verify_bundle()
    state=a.launch_state.resolve();state.mkdir(parents=True,exist_ok=True)
    identity=dict(project=str(a.project.resolve()),storage_state=str(a.storage_state.resolve()),
                  branch=str(a.branch.resolve()),bundle=str(BUNDLE),target=40)
    with (state/'.launch.lock').open('a') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        ident=state/'identity.json'
        if ident.exists() and json.loads(ident.read_text())!=identity: raise RuntimeError('Launch state belongs to another run')
        storage.atomic_json(ident,identity)
        current=os.environ.get('SLURM_JOB_ID','') if a.mode=='gate' else ''
        rows=storage.queue(); waiting=classify(rows,current)
        if a.mode=='submit':
            if (state/'last_submission.json').exists() and not a.resume:
                raise RuntimeError('Launch receipt already exists. Use resume only after checking the logs and queue')
            if len(rows)>=300:
                raise RuntimeError('All 300 user job slots are occupied. Retry this command after a slot becomes free; nothing was submitted')
            submit_job(a,'gate',delay=True); return
        if len([r for r in rows if r[0]==current and r[1]=='slide_roll_gate'])!=1:
            raise RuntimeError('Gate must run as its submitted Slurm job')
        if waiting:
            print('Existing slider jobs still active; leaving them unchanged: '+repr(waiting[:8]),flush=True)
            submit_job(a,'gate',dependency=current,delay=True); return
        compact=json.loads((a.storage_state.resolve()/'COMPLETE.json').read_text())
        if compact.get('project')!=str(a.project.resolve()) or 'saved_bytes' not in compact:
            raise RuntimeError('Storage compaction has not completed')
        if not (a.branch.resolve()/'INITIALIZED.json').exists():
            with storage.maintenance_lock(a.project.resolve()):
                storage.check_sources(a.project.resolve())
                subprocess.check_call([os.environ.get('RSCRIPT','Rscript'),'--vanilla',str(BUNDLE/'check_source.R'),str(BUNDLE),str(a.project.resolve())])
        submit_job(a,'setup',dependency=current)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['submit','gate'])
    for name in ('project','storage_state','branch','launch_state'): p.add_argument(name,type=Path)
    p.add_argument('--resume',action='store_true');a=p.parse_args()
    try: run(a)
    except Exception as e:
        print('STOP: '+str(e)+'. Existing jobs were not cancelled.',file=sys.stderr);sys.exit(1)
