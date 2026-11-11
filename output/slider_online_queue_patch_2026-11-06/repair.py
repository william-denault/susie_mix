"""Install retryable scheduler reads and resume the validated batch-6 commit."""
import hashlib, os, shutil, subprocess, sys, time
from pathlib import Path
import repair_support as support
from repair_support import read,write,encoded,atomic,digest,inside,move
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
PROJECT=Path('/project2/mstephens/wdenault/susie_mix')
TAG='scheduler_queue_retry_2026-11-06'

def confirm_failure(branch,project):
    d=project/'results_slide_prior_em/iteration_015/batches/batch_06'
    error=read(branch/'ERROR.json')
    if error.get('job')!='49863377':raise RuntimeError('Different stopping controller; inspect before repair')
    command=['squeue','--noheader','--array','--states=all','--user','wdenault','--format=%i|%F|%100j|%T|%r|%S']
    if error.get('error')!=str(subprocess.TimeoutExpired(command,15)):
        raise RuntimeError('Not the reported squeue timeout; refusing to clear another error')
    if (d/'commit.rds').exists():raise RuntimeError('Batch 6 already has a commit; refusing a stage-5 repair')
    return d

def repair(project=PROJECT,expected_batches=(2030,2030,2029,2029,2029,2029)):
    project=project.resolve();branch=project/'workflows/slider_online_015_025';inside(branch,project)
    support.HERE=HERE;mb=support.runtime()
    support.verify_manifest(HERE,(HERE/'SHA256SUMS').read_bytes())
    old=(HERE/'baseline_SHA256SUMS').read_bytes();new=(HERE/'runtime/SHA256SUMS').read_bytes()
    state=branch/'repairs'/TAG;state.mkdir(parents=True,exist_ok=True)
    with mb.lock(branch/'.queue_repair.lock'):
      if (state/'COMPLETE.json').exists():
        print('Scheduler-query repair already complete; checking continuation.',flush=True)
      else:
       with mb.lock(branch/'.launch.lock'),mb.lock(branch/'.controller.lock'):
        support.idle(mb)
        if (branch/'pending_submission.json').exists():raise RuntimeError('Unresolved submission; refusing to alter it')
        planfile=state/'PLAN.json'
        if not planfile.exists():
            print('Checking five committed updates and all 12176 saved results.',flush=True)
            support.verify_manifest(branch/'code',old);mb.verify_inputs(branch)
            cfg=read(branch/'config.json')
            expected={'schema':3,'source_iteration':14,'rho':.3,'passes':6,'batches':9,'full_steps':5,'workers':298,'poll_minutes':1}
            if any(cfg.get(k)!=v for k,v in expected.items()) or Path(cfg['source']).resolve()!=project or \
               cfg.get('bundle_sha256')!=hashlib.sha256(old).hexdigest():raise RuntimeError('Unexpected installed workflow')
            if int(mb.csvrows(branch/'progress.csv')[0]['stage'])!=5:raise RuntimeError('Expected five committed stages')
            confirm_failure(branch,project)
            protected=support.protected_state(branch,mb)
            if tuple(protected['validated_per_batch'])!=expected_batches:raise RuntimeError('Unexpected saved batch sizes')
            before=state/'before';before.mkdir(exist_ok=True)
            for rel in ('config.json','frozen_inputs.json','ERROR.json'):
                shutil.copyfile(str(branch/rel),str(before/rel))
            cfg['bundle_sha256']=hashlib.sha256(new).hexdigest();cfgbytes=encoded(cfg)
            frozen=read(branch/'frozen_inputs.json');frozen['config.json']=hashlib.sha256(cfgbytes).hexdigest()
            write(state/'config.after.json',cfgbytes);atomic(state/'frozen.after.json',frozen)
            atomic(planfile,dict(project=str(project),branch=str(branch),new_manifest=hashlib.sha256(new).hexdigest(),
                protected=protected,expected_batches=list(expected_batches),time=time.time()))
        plan=read(planfile)
        if plan['project']!=str(project) or plan['branch']!=str(branch) or \
           plan['new_manifest']!=hashlib.sha256(new).hexdigest():raise RuntimeError('Foreign repair transaction')
        support.check_protected(plan,mb)
        if not (state/'PREFLIGHT_OK.json').exists():
            print('Dry-running the batch-6 commit with validated summaries; no fitting.',flush=True)
            subprocess.check_call([os.environ.get('RSCRIPT','Rscript'),'--vanilla',str(HERE/'check_commit.R'),
                str(HERE/'runtime'),str(branch),str(state),str(plan['expected_batches'][5])])
            support.check_protected(plan,mb)
            atomic(state/'PREFLIGHT_OK.json',dict(input_state_unchanged=True,dry_run_stage=6,time=time.time()))
        support.idle(mb);support.patch_code(branch,state,old,new)
        for rel,after in [('config.json','config.after.json'),('frozen_inputs.json','frozen.after.json')]:
            original=(state/'before'/rel).read_bytes();updated=(state/after).read_bytes()
            if (branch/rel).read_bytes() not in (original,updated):raise RuntimeError('Changed configuration during repair')
            write(branch/rel,updated)
        error=branch/'ERROR.json';archive=state/'ERROR.before_repair.json';expected=(state/'before/ERROR.json').read_bytes()
        if error.exists():
            if error.read_bytes()!=expected:raise RuntimeError('A different failure appeared during repair')
            move(error,archive,branch)
        elif not archive.exists() or archive.read_bytes()!=expected:raise RuntimeError('Missing error backup')
        mb.HASH_CACHE.clear();mb.verify_code(branch/'code');mb.verify_inputs(branch);support.check_protected(plan,mb)
        atomic(state/'COMPLETE.json',dict(preserved_validated_results=sum(plan['expected_batches']),
            committed_stages_before_resume=5,resume_iteration=15,resume_batch=6,
            dry_run_commit_passed=True,existing_results_and_priors_unchanged=True,time=time.time()))
        print('Scheduler-query repair verified: %d validated results preserved; five prior updates unchanged.'%sum(plan['expected_batches']),flush=True)
        print('Resuming the batch-6 commit, then batch 7; no completed gene needs refitting.',flush=True)
      support.resume(branch,mb)

if __name__=='__main__':
    try:repair()
    except Exception as exc:print('STOP: '+str(exc),file=sys.stderr);sys.exit(1)
