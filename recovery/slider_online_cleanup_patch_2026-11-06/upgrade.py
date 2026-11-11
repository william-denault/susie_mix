"""Pause one pending controller, drain its workers, upgrade code, release it."""
import contextlib, hashlib, importlib.util, os, re, shutil, subprocess, sys, time
from pathlib import Path
import support
from support import read, atomic, encoded, write, digest, move
HERE=Path(__file__).resolve().parent
BRANCH=Path('/project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025')
TAG='cleanup_transition_2026-11-06'

def runtime():
    support.HERE=HERE
    return support.runtime()

def state_dir(branch):return branch/'repairs'/TAG

def details(job):
    mb=runtime();command=['scontrol','show','job','-o',str(job)]
    if os.environ.get('SLURM_JOB_ID'):
        with mb.wait_for_scheduler():text=mb.queue_output(command,os.environ.copy())
    else:text=mb.queue_output(command,os.environ.copy())
    return dict(re.findall(r'(\w+)=(\S+)',text))

def control_action(action,job):
    if action not in ('hold','release'):raise RuntimeError('Unsupported controller action')
    return subprocess.run(['scontrol',action,job],check=True,stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,universal_newlines=True,timeout=60)

def identity(info,branch,mb):
    cfg=read(branch/'config.json')
    if info.get('JobName')!=mb.names(cfg)[0] or info.get('UserId','').split('(')[0]!=os.environ['USER'] or \
       info.get('WorkDir')!=str(branch) or info.get('Command')!=str(branch/'code/controller_job'):
        raise RuntimeError('Controller identity does not match this workflow')

def release(branch,mb):
    state=state_dir(branch);target=state/'hold.json'
    if not target.exists():return
    held=read(target);job=held['job']
    while True:
        rows=mb.queue()
        if not any(r['job']==job for r in rows):
            atomic(state/'release.json',dict(job=job,already_absent=True,time=time.time()));return
        try:info=details(job)
        except subprocess.CalledProcessError:
            if not any(r['job']==job for r in mb.queue()):continue
            raise
        identity(info,branch,mb)
        if info.get('Priority')!='0':break
        try:control_action('release',job)
        except (subprocess.TimeoutExpired,subprocess.CalledProcessError,OSError) as exc:
            if not mb.transient_queue_error(exc):raise
            print('Release response delayed; checking the exact controller before retrying.',flush=True)
            time.sleep(20)
        # Read back state; a lost acknowledgement must not cause a blind action.
    atomic(state/'release.json',dict(job=job,time=time.time()))
    print('Released existing controller '+job+'; completed genes and queued work are reused.',flush=True)

def validate_baseline(branch,mb):
    support.verify_manifest(branch/'code',(HERE/'baseline_SHA256SUMS').read_bytes())
    mb.verify_inputs(branch)
    cfg=read(branch/'config.json')
    if cfg.get('bundle_sha256')!=digest(HERE/'baseline_SHA256SUMS') or \
       any(cfg.get(k)!=v for k,v in dict(source_iteration=14,rho=.3,passes=6,batches=9,full_steps=5).items()):
        raise RuntimeError('Unexpected installed code/configuration; no change made')
    if (branch/'ERROR.json').exists() or (branch/'pending_submission.json').exists():
        raise RuntimeError('Existing failure or ambiguous submission: inspect before upgrading')

def protected(branch,mb):
    root=Path(read(branch/'config.json')['source'])/'results_slide_prior_em'
    paths={branch/n for n in ('state.rds','model.rds','stages.csv','progress.csv','prior_history.csv','latest_priors.csv','original_prior_history.csv')}
    paths.update((branch/'plan').glob('*.txt'))
    paths.update(root/n for n in ('prior_history.csv','prior_history_long.csv','online_latest_priors.csv'))
    for d in root.glob('iteration_0*'):
        if d.name>='iteration_015':
            paths.update(d.glob('*.csv'))
            for b in (d/'batches').glob('*'):
                paths.update(b.glob('*.csv'));paths.update(b.glob('*.rds'))
    return {str(p):digest(p) for p in paths if p.is_file()}

def verify_protected(plan):
    for p,wanted in plan['protected'].items():
        if digest(p)!=wanted:raise RuntimeError('Analysis changed while paused: '+p)

def drained(branch,mb):
    rows=mb.queue();own=mb.names(read(branch/'config.json'));job=read(state_dir(branch)/'hold.json')['job']
    if mb.legacy_jobs(rows):raise RuntimeError('Legacy slider jobs reappeared')
    target=[r for r in rows if r['job']==job]
    if len(target)!=1 or target[0]['name']!=own[0] or target[0]['state']!='PENDING':
        raise RuntimeError('Held controller changed; refusing to install')
    info=details(job);identity(info,branch,mb)
    if info.get('Priority')!='0':raise RuntimeError('Controller is no longer held')
    return not any(r['name'] in own and r['job']!=job for r in rows)

def install(branch,mb):
    state=state_dir(branch);planfile=state/'PLAN.json'
    old=(HERE/'baseline_SHA256SUMS').read_bytes();new=(HERE/'runtime/SHA256SUMS').read_bytes()
    if not drained(branch,mb):raise RuntimeError('Workers/controllers still active; no code swap')
    with mb.lock(branch/'.launch.lock'),mb.lock(branch/'.controller.lock'):
        if not planfile.exists():
            validate_baseline(branch,mb)
            before=state/'before';before.mkdir(exist_ok=True)
            for name in ('config.json','frozen_inputs.json'):
                shutil.copyfile(str(branch/name),str(before/name))
            cfg=read(branch/'config.json');cfg['bundle_sha256']=hashlib.sha256(new).hexdigest()
            config_bytes=encoded(cfg);write(state/'config.after.json',config_bytes)
            frozen=read(branch/'frozen_inputs.json');frozen['config.json']=hashlib.sha256(config_bytes).hexdigest()
            atomic(state/'frozen.after.json',frozen)
            atomic(planfile,dict(branch=str(branch),new_manifest=hashlib.sha256(new).hexdigest(),
                protected=protected(branch,mb),time=time.time()))
        plan=read(planfile)
        if plan['branch']!=str(branch) or plan['new_manifest']!=hashlib.sha256(new).hexdigest():
            raise RuntimeError('Different upgrade transaction')
        verify_protected(plan)
        support.patch_code(branch,state,old,new)
        for rel,after in [('config.json','config.after.json'),('frozen_inputs.json','frozen.after.json')]:
            original=(state/'before'/rel).read_bytes();updated=(state/after).read_bytes()
            if (branch/rel).read_bytes() not in (original,updated):raise RuntimeError('Configuration changed during upgrade')
            write(branch/rel,updated)
        mb.HASH_CACHE.clear();mb.verify_code(branch/'code');mb.verify_inputs(branch);verify_protected(plan)
        atomic(state/'INSTALLED.json',dict(time=time.time(),analysis_preserved=True,worker_jobs_cancelled=0))
        print('Cleanup hardening installed; all fitted results, batch updates and history preserved.',flush=True)

def recover(branch,mb):
    state=state_dir(branch)
    if (state/'PLAN.json').exists() and not (state/'INSTALLED.json').exists():install(branch,mb)
    # With no PLAN, installed code was never touched; release the pending job.
    release(branch,mb)

def gate(branch=BRANCH):
    mb=runtime();state=state_dir(branch);state.mkdir(parents=True,exist_ok=True)
    support.verify_manifest(HERE,(HERE/'SHA256SUMS').read_bytes())
    with mb.wait_for_scheduler(),mb.lock(state/'gate.lock'):
        if (state/'INSTALLED.json').exists():release(branch,mb);return
        if (state/'hold.json').exists():
            # A restarted gate first restores any interrupted transaction.
            recover(branch,mb)
            if not (state/'INSTALLED.json').exists():raise RuntimeError('Previous attempt released; inspect its log before retrying')
            return
        validate_baseline(branch,mb)
        launch=read(HERE/'launch.json')
        if launch.get('gate')!=os.environ.get('SLURM_JOB_ID') or not launch.get('watchdog'):
            raise RuntimeError('No confirmed recovery job; refusing to hold a controller')
        watchdog=details(launch['watchdog'])
        if watchdog.get('JobName')!='slide_cleanup_recover' or watchdog.get('UserId','').split('(')[0]!=os.environ['USER'] or \
           watchdog.get('WorkDir')!=str(HERE) or watchdog.get('Command')!=str(HERE/'recover_job') or \
           ('afterany:'+launch['gate']) not in watchdog.get('Dependency',''):
            raise RuntimeError('Recovery-job identity/dependency mismatch')
        deadline=time.monotonic()+18*3600
        try:
            while not (state/'hold.json').exists():
                if (branch/'DONE.json').exists():print('Analysis already complete; no upgrade required.',flush=True);return
                rows=mb.queue();ctl=mb.names(read(branch/'config.json'))[0]
                candidates=[r for r in rows if r['name']==ctl and r['state']=='PENDING']
                if not candidates and not any(r['name']==ctl for r in rows):
                    raise RuntimeError('No live controller; inspect workflow status')
                if candidates:
                    if len(candidates)!=1:raise RuntimeError('Multiple pending controllers; refusing to guess')
                    job=candidates[0]['job'];info=details(job);identity(info,branch,mb)
                    if info.get('JobState')=='PENDING' and info.get('Priority')!='0':
                        atomic(state/'hold.json',dict(job=job,time=time.time(),original_priority=info['Priority']))
                        # If the RPC fails after being accepted, recovery still knows the exact job.
                        control_action('hold',job)
                        info=details(job);identity(info,branch,mb)
                        if info.get('JobState')!='PENDING' or info.get('Priority')!='0':
                            raise RuntimeError('Controller started before hold; leaving its work unchanged')
                        print('Temporarily held pending controller '+job+'; every gene job continues.',flush=True)
                        break
                    if info.get('Priority')=='0':raise RuntimeError('Controller was already held by another action')
                if time.monotonic()>deadline:raise RuntimeError('No pending-controller window within 18 hours')
                time.sleep(20)
            while not drained(branch,mb):
                if time.monotonic()>deadline:raise RuntimeError('Drain deadline reached; restoring normal continuation')
                print('Waiting for already-submitted gene jobs/current controller to finish; no jobs cancelled.',flush=True)
                time.sleep(30)
            install(branch,mb)
        finally:
            # The afterany recovery job also executes this after TIMEOUT/node failure.
            recover(branch,mb)

def launch():
    mb=runtime();support.verify_manifest(HERE,(HERE/'SHA256SUMS').read_bytes())
    ledger=HERE/'launch.json';pending=HERE/'launcher_pending.json'
    with mb.lock(HERE/'.launch.lock'):
        if ledger.exists() or pending.exists():raise RuntimeError('Upgrade was already submitted or has an uncertain receipt; inspect launch.json/logs, do not duplicate it')
        if len(mb.queue())>298:raise RuntimeError('Two free submission slots are needed; current EM is unchanged. Retry after some jobs finish.')
        (HERE/'logs').mkdir(exist_ok=True)
        cmd=['sbatch','--parsable','--exclude=midway2-0003,midway2-0011','--hold',
             '--chdir='+str(HERE),'--output='+str(HERE/'logs/gate_%j.out'),'--error='+str(HERE/'logs/gate_%j.err'),
             str(HERE/'gate_job'),str(HERE)]
        atomic(pending,dict(kind='gate',command=cmd))
        job=subprocess.check_output(cmd,universal_newlines=True,timeout=120).strip().split(';')[0]
        if not job.isdigit():raise RuntimeError('Ambiguous gate submission')
        atomic(ledger,dict(gate=job));pending.unlink()
        cmd=['sbatch','--parsable','--exclude=midway2-0003,midway2-0011','--dependency=afterany:'+job,
             '--chdir='+str(HERE),'--output='+str(HERE/'logs/recover_%j.out'),'--error='+str(HERE/'logs/recover_%j.err'),
             str(HERE/'recover_job'),str(HERE)]
        atomic(pending,dict(kind='watchdog',command=cmd))
        watchdog=subprocess.check_output(cmd,universal_newlines=True,timeout=120).strip().split(';')[0]
        if not watchdog.isdigit():raise RuntimeError('Ambiguous recovery-job submission')
        atomic(ledger,dict(gate=job,watchdog=watchdog));pending.unlink()
        control_action('release',job)
        print('Submitted upgrade gate '+job+' and recovery job '+watchdog+'. No gene jobs are cancelled.',flush=True)

if __name__=='__main__':
    try:
        if sys.argv[1]=='launch':launch()
        elif sys.argv[1]=='gate':gate()
        elif sys.argv[1]=='recover':
            mb=runtime();state=state_dir(BRANCH);state.mkdir(parents=True,exist_ok=True)
            with mb.wait_for_scheduler(),mb.lock(state/'gate.lock'):recover(BRANCH,mb)
        else:raise RuntimeError('Unknown mode')
    except Exception as exc:print('STOP: '+str(exc),file=sys.stderr);sys.exit(1)
