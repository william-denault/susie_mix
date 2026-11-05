"""In-place handoff after 014: online EM 015..020, full EM 021..025."""
import argparse, contextlib, csv, datetime, hashlib, json, os, random, re
import shutil, subprocess, sys, time
from pathlib import Path
sys.dont_write_bytecode = True
SOURCE = '/project2/mstephens/wdenault/susie_mix'
BRANCH = SOURCE+'/workflows/slider_online_015_025'
TERMINAL = {'COMPLETED','CANCELLED','FAILED','TIMEOUT','NODE_FAIL','OUT_OF_MEMORY',
            'PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED'}
HASH_CACHE = {}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def atomic(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp.'+str(os.getpid()))
    with tmp.open('w') as f:
        json.dump(obj,f,indent=2,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(str(tmp),str(path))

def read(path): return json.loads(Path(path).read_text())
def cached_sha(path):
    path=Path(path);s=path.stat();key=(str(path),s.st_size,s.st_mtime_ns)
    if key not in HASH_CACHE: HASH_CACHE[key]=sha(path)
    return HASH_CACHE[key]
def csvrows(path):
    with Path(path).open(newline='') as f: return list(csv.DictReader(f))
def writecsv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    with tmp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    os.replace(str(tmp),str(path))

@contextlib.contextmanager
def lock(path):
    import fcntl
    with Path(path).open('a') as f:
        fcntl.flock(f.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        yield

def verify_code(code):
    code=Path(code)
    for line in (code/'SHA256SUMS').read_text().splitlines():
        digest,name=line.split('  ',1)
        if sha(code/name)!=digest: raise RuntimeError('Bundle file changed: '+name)

def verify_source(source):
    rt=source/'slider_inplace_runtime_2026-10-29'
    frozen=read(rt/'frozen.json')
    if Path(frozen['project']).resolve()!=source: raise RuntimeError('Source runtime belongs to another project')
    for rel,digest in frozen['hashes'].items():
        p=(source/rel).resolve()
        if source not in p.parents or sha(p)!=digest: raise RuntimeError('Source runtime changed: '+rel)

def rrun(branch,mode,*args):
    code=branch/'code'
    cmd=[os.environ.get('RSCRIPT','Rscript'),'--vanilla',str(code/'engine.R'),mode,str(code)]
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    subprocess.check_call(cmd+[str(a) for a in args],env=env)

def queue():
    env=os.environ.copy();env['TZ']='UTC'
    out=subprocess.check_output(['squeue','--noheader','--array','--states=all',
        '--user',os.environ['USER'],'--format=%i|%F|%100j|%T|%r|%S'],env=env,universal_newlines=True,timeout=15)
    rows=[]
    for line in out.splitlines():
        if not line.strip(): continue
        fields=[s.strip() for s in line.split('|')]
        if len(fields)!=6: raise RuntimeError('Invalid queue response')
        job,array,name,state,reason,start=fields
        # %r is the raw reason. Also accept %R-style parentheses from older
        # fixtures/wrappers: comparing '(BeginTime)' to 'BeginTime' blocked
        # the original handoff indefinitely.
        if reason.startswith('(') and reason.endswith(')'):reason=reason[1:-1]
        if state not in TERMINAL: rows.append(dict(job=job,array=array,name=name,state=state,reason=reason,start=start))
    if len(set(r['job'] for r in rows))!=len(rows): raise RuntimeError('Duplicate queue entries')
    return rows

def names(cfg): return 'slide_online_c_'+cfg['id'], 'slide_online_f_'+cfg['id']
def legacy_jobs(rows):
    return [r for r in rows if r['name'].startswith(('slide_inplace_','slide_parallel_')) or
            r['name'] in ('slide_next_gate','slide_next_setup','em_susie_slide_prior','em_slide_prior_prepare')]

def result_root(branch): return Path(read(branch/'config.json')['source'])/'results_slide_prior_em'

def assert_owned_iteration(branch,iteration):
    root=result_root(branch);d=root/('iteration_%03d'%iteration)
    cfg=read(branch/'config.json');owner=d/'ONLINE_OWNER.json'
    if d.is_symlink() or d.resolve().parent!=root.resolve():raise RuntimeError('Invalid iteration path')
    if d.exists() and (not owner.exists() or read(owner).get('workflow_id')!=cfg['id']):
        raise RuntimeError('Iteration already exists outside this workflow: '+str(d))
    if not d.exists():
        staging=root/('.online_prepare_%03d_%s'%(iteration,cfg['id']))
        staging.mkdir(exist_ok=True)
        atomic(staging/'ONLINE_OWNER.json',dict(workflow_id=cfg['id'],iteration=iteration,source_iteration=14))
        os.rename(str(staging),str(d))
def external_submitters(rows,cfg):
    own=names(cfg)
    return [r for r in rows if r['name'] not in own and
            (r['name'].endswith('control') or 'launch' in r['name'] or
             r['name'] in ('slide_next_gate','slide_next_setup','em_susie_mix','em_slide_prior_prepare'))]

def submission_budget(rows,cfg):
    # The original controller exits if it cannot queue its successor. Keep one
    # additional slot per external controller; its queued job already counts
    # in rows, but its not-yet-submitted successor does not.
    return 300-len(external_submitters(rows,cfg))

def safe_window(rows,cfg,now=None):
    """Never take the original controller's nonblocking .prepare.lock.

    Its queue snapshot must not be in flight while we reserve slots. A future
    BeginTime is a hard earliest start, unlike a Priority start-time estimate.
    All our sbatch RPCs have a 10-second deadline within this 90-second margin.
    """
    now=time.time() if now is None else now
    for r in external_submitters(rows,cfg):
        if r['state']!='PENDING' or r['reason']!='BeginTime': return False
        try:
            earliest=datetime.datetime.strptime(r['start'],'%Y-%m-%dT%H:%M:%S').replace(tzinfo=datetime.timezone.utc).timestamp()
        except ValueError: return False
        if earliest-now<90: return False
    return True

def submission_queue(branch,needed=1,wait_seconds=600):
    cfg=read(branch/'config.json');deadline=time.monotonic()+wait_seconds
    while True:
        rows=queue()
        if safe_window(rows,cfg) and len(rows)+needed<=submission_budget(rows,cfg): return rows
        if time.monotonic()>=deadline: raise RuntimeError('No safe submission window/free slots; retry --resume later')
        print('Waiting for free slots and an idle original controller; existing jobs remain unchanged.',flush=True)
        time.sleep(20)

def submit(branch,kind,rows,task=None,dependency=None):
    cfg=read(branch/'config.json');ctl,worker=names(cfg)
    pending=branch/'pending_submission.json'
    if pending.exists(): raise RuntimeError('Ambiguous previous submission; inspect '+str(pending))
    if not safe_window(rows,cfg): raise RuntimeError('Unsafe submission window')
    count=len(task['genes']) if task else 1
    if len(rows)+count>submission_budget(rows,cfg):
        raise RuntimeError('Submission would consume the 300-job budget or an external controller successor slot')
    logs=Path(task['directory']) if task else branch/'logs';logs.mkdir(parents=True,exist_ok=True)
    cmd=['sbatch','--parsable','--export=ALL','--time=23:00:00',
         '--exclude=midway2-0003,midway2-0011','--chdir='+str(branch),
         '--job-name='+(worker if task else ctl),
         '--output='+str(logs/('worker_%A_%a.out' if task else 'control_%j.out')),
         '--error='+str(logs/('worker_%A_%a.err' if task else 'control_%j.err'))]
    if task: cmd+=['--array=1-'+str(count)]
    else:
        cmd+=['--begin=now+1minutes']
        if dependency: cmd+=['--dependency=afterany:'+dependency]
    cmd += [str(branch/'code'/('worker_job' if task else 'controller_job')),str(branch)]
    if task: cmd.append(str(Path(task['directory'])/'task.json'))
    payload=dict(kind=kind,command=cmd,time=time.time(),task=task)
    atomic(pending,payload)
    # An interrupted/ambiguous RPC is deliberately never retried automatically.
    job=subprocess.check_output(cmd,universal_newlines=True,timeout=10).strip().split(';')[0]
    if not job.isdigit(): raise RuntimeError('Ambiguous sbatch response; inspect pending_submission.json')
    payload['job']=job
    if task: atomic(Path(task['directory'])/'submission.json',payload)
    else: atomic(branch/'last_controller.json',payload)
    pending.unlink()
    print('Submitted '+kind+' '+job,flush=True)
    return job

def balanced(genes,groups,seed):
    if len(genes)<groups or len(set(genes))!=len(genes): raise ValueError('Invalid balanced cohort')
    shuffled=list(genes);random.Random(seed).shuffle(shuffled)
    q,r=divmod(len(genes),groups);out=[];pos=0
    for b in range(groups):
        n=q+(b<r);out.append(shuffled[pos:pos+n]);pos+=n
    return out

def make_plan(branch,genes):
    cfg=read(branch/'config.json');stages=[];sid=0
    for epoch in range(1,7):
        batches=balanced(genes,9,cfg['seed']+epoch)
        for batch,selected in enumerate(batches,1):
            sid+=1;iteration=14+epoch;group='../../results_slide_prior_em/iteration_%03d'%iteration
            directory=group+'/batches/batch_%02d'%batch
            gf='plan/pass_%03d_batch_%02d.txt'%(epoch,batch)
            p=branch/gf;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('\n'.join(selected)+'\n')
            stages.append(dict(stage=sid,phase='minibatch',pass_=epoch,batch=batch,rho=0.3,
                               directory=directory,genes_file=gf,group=group,iteration=iteration,results_dir=group+'/results'))
    for step in range(1,6):
        sid+=1;gf='plan/full_%03d.txt'%step;(branch/gf).write_text('\n'.join(genes)+'\n')
        iteration=20+step;group='../../results_slide_prior_em/iteration_%03d'%iteration
        stages.append(dict(stage=sid,phase='full',pass_=step,batch=1,rho=1.,
             directory=group+'/batches/full_step',genes_file=gf,group=group,iteration=iteration,results_dir=group+'/results'))
    stages=[{('pass' if k=='pass_' else k):v for k,v in row.items()} for row in stages]
    writecsv(branch/'stages.csv',stages)
    return stages

def pin(expected,results):
    """Hard links, not symlinks: original cleanup cannot invalidate the seed."""
    results.mkdir(parents=True,exist_ok=True)
    for r in expected:
        if not re.fullmatch(r'[A-Za-z0-9_.-]+',r['gene']) or r['gene'] in ('.','..'): raise RuntimeError('Invalid gene name')
        src=Path(r['source']);dst=results/(r['gene']+'.rds')
        if not dst.exists(): os.link(str(src),str(dst),follow_symlinks=False)
        st=dst.stat()
        if dst.is_symlink() or st.st_size!=int(float(r['size'])) or abs(st.st_mtime-float(r['mtime']))>0.001:
            raise RuntimeError('Source changed while pinning: '+r['gene'])

def release_results(branch,directory):
    directory=Path(directory)
    root=branch.resolve();d=directory.resolve()
    if root not in d.parents or directory.is_symlink() or directory.name!='results':
        raise RuntimeError('Invalid branch-only cleanup target: '+str(directory))
    if not directory.exists(): return
    files=list(directory.iterdir())
    if any(p.is_symlink() or not p.is_file() or p.suffix!='.rds' for p in files):
        raise RuntimeError('Unexpected object in branch results; nothing removed')
    for p in files: p.unlink()

def seed_run(branch):
    cfg=read(branch/'config.json');source=Path(cfg['source'])
    ready=branch/'seed/READY.json'
    if ready.exists():
        if legacy_jobs(queue()):raise RuntimeError('Legacy slider jobs reappeared after checkpoint capture')
        snapshot=Path(read(ready)['snapshot'])
    else:
        if legacy_jobs(queue()):
            print('Waiting for the existing jobs to finish and certify iteration 014.',flush=True)
            return False
        if not (source/'results_slide_prior_em/iteration_014/recovery/parallel/COMPLETE.rds').exists():
            raise RuntimeError('No legacy jobs remain, but 014 is not certified. Inspect finisher logs before resuming.')
        verify_source(source)
        root=branch/'seed';root.mkdir(exist_ok=True)
        attempts=sorted(root.glob('snapshot_*'))
        snapshot=attempts[-1] if attempts else root/'snapshot_001'
        if (snapshot/'CAPTURE_FAILED.json').exists():
            snapshot=root/('snapshot_%03d'%(len(attempts)+1))
        # A partially pinned checkpoint may be reclaimed by original cleanup.
        # Retry capture from its latest certified successor, without touching it.
        for attempt in range(3):
            try:
                if not (snapshot/'capture.csv').exists(): rrun(branch,'capture',source,snapshot)
                expected=csvrows(snapshot/'expected.csv')
                pin(expected,snapshot/'results')
                capture=csvrows(snapshot/'capture.csv')[0]
                if int(capture['iteration'])!=14:raise RuntimeError('Handoff requires exactly completed iteration 014')
                atomic(ready,dict(snapshot=str(snapshot),source=str(source),time=time.time(),
                     source_iteration=int(capture['iteration']),genes=int(capture['genes']),
                     tissue_fits=int(capture['tissue_fits'])))
                break
            except (FileNotFoundError,subprocess.CalledProcessError) as exc:
                atomic(snapshot/'CAPTURE_FAILED.json',dict(error=str(exc),time=time.time()))
                if isinstance(exc,subprocess.CalledProcessError):
                    if exc.returncode!=75: raise
                    print('Source certificate is being refreshed; the queued successor will retry capture.',flush=True)
                    return False
                if attempt==2:
                    print('Source cleanup raced checkpoint capture; the queued successor will retry.',flush=True)
                    return False
                snapshot=root/('snapshot_%03d'%(len(list(root.glob('snapshot_*')))+1))
        for old in root.glob('snapshot_*'):
            if old!=snapshot: release_results(branch,old/'results')
    if not (branch/'stages.csv').exists(): make_plan(branch,[r['gene'] for r in csvrows(snapshot/'genes.csv')])
    rrun(branch,'initialize',branch,snapshot)
    frozen=branch/'frozen_inputs.json'
    if not frozen.exists():
        files=[branch/'config.json',branch/'model.rds',branch/'stages.csv',branch/'original_prior_history.csv',
               snapshot/'source_priors.csv',snapshot/'metadata.rds',snapshot/'expected.csv']
        files+=list((branch/'plan').glob('*.txt'))
        atomic(frozen,{p.relative_to(branch).as_posix():sha(p) for p in files})
    verify_inputs(branch)
    rrun(branch,'history',branch)
    return True

def finish_source(branch):
    """Finish 014 through the unchanged private planner and gene workers.

    One controller now owns refills, validation and the online continuation.
    Existing workers keep their original Slurm IDs and completion receipts.
    """
    if (branch/'seed/READY.json').exists():return True
    cfg=read(branch/'config.json');source=Path(cfg['source']);root=source/'results_slide_prior_em'
    if not (branch/'HANDOFF.json').exists():raise RuntimeError('Source finish requires a committed handoff')
    import handoff
    handoff.boundary(source,sys.modules[__name__])
    verify_source(source)
    rows=submission_queue(branch,0,20*3600)
    handoff.check_queue(rows,cfg,sys.modules[__name__])
    old=legacy_jobs(rows)
    if any(r['name'] not in ('slide_inplace_fit','slide_inplace_validate') for r in old):
        raise RuntimeError('A legacy controller reappeared; refusing competing planning')
    capacity=max(0,min(cfg['workers']-len(old),submission_budget(rows,cfg)-len(rows)))
    rt=source/'slider_inplace_runtime_2026-10-29'
    queue_file=branch/'logs'/('finish_014_queue_'+os.environ['SLURM_JOB_ID']+'.txt')
    queue_file.write_text(''.join('%s|%s|%s\n'%(r['job'],r['array'],r['name']) for r in rows))
    env=os.environ.copy();env.update(SUSIE_MIX_PROJECT_DIR=str(source),SUSIE_SLIDE_TARGET_ITERATION='14',
        SUSIE_SLIDE_WORKERS=str(cfg['workers']),SUSIE_SLIDE_POLL_MINUTES='1',
        SUSIE_SLIDE_WAIT_FOR_SIMULATIONS='0',SUSIE_SLIDE_MAX_ITER='5000',SUSIE_SLIDE_RECORD_ELBO='0',
        SUSIE_SLIDE_EXCLUDE_NODES='midway2-0003,midway2-0011',SBATCH_TIMELIMIT='23:00:00')
    with lock(root/'.prepare.lock'):
        if (root/'parallel_pending_submission.txt').exists():raise RuntimeError('Unresolved legacy submission; inspect before continuing')
        command=[os.environ.get('RSCRIPT','Rscript'),'--vanilla',
                 str(rt/'script/scan_tissue_attempt/run_slide_prior_parallel.R'),
                 'plan',str(source),'14',str(capacity),str(queue_file)]
        response=subprocess.check_output(command,env=env,universal_newlines=True).strip().splitlines()
        if len(response)<2 or response[1]!='14':raise RuntimeError('Unexpected source-014 plan response')
        status=response[0]
        if status=='complete':
            if old:raise RuntimeError('014 marked complete while legacy workers still active')
            if not (root/'iteration_014/recovery/parallel/COMPLETE.rds').exists():raise RuntimeError('014 completion certificate missing')
            atomic(branch/'SOURCE_014_COMPLETE.json',dict(iteration=14,time=time.time()))
            print('Original iteration 014 certified; continuing to online 015.',flush=True)
            return True
        if status=='waiting':
            print('Finishing original 014: %d active fitting/validation tasks; waiting for completions.'%len(old),flush=True)
            return False
        if status!='batch' or len(response)!=5 or not response[3].isdigit() or response[4] not in ('fit','validate'):
            raise RuntimeError('Invalid source-014 worker plan')
        batch=Path(response[2]).resolve();count=int(response[3]);action=response[4]
        expected=(root/'iteration_014/recovery/parallel/batches').resolve()
        if batch.parent!=expected or not batch.is_dir() or not 1<=count<=capacity:
            raise RuntimeError('Source worker plan is outside the expected directory/capacity')
        pending=branch/'pending_submission.json'
        if pending.exists():raise RuntimeError('Unresolved submission; inspect pending_submission.json')
        command=['sbatch','--parsable','--export=ALL','--time=23:00:00',
                 '--exclude=midway2-0003,midway2-0011','--array=1-'+str(count),
                 '--job-name=slide_inplace_'+action,'--chdir='+str(rt/'job'),
                 '--output='+str(batch/'worker_%A_%a.out'),'--error='+str(batch/'worker_%A_%a.err'),
                 str(rt/'job/slide_prior_parallel_worker'),str(source),'14',str(batch)]
        payload=dict(kind='finish_014_workers',command=command,directory=str(batch),time=time.time())
        atomic(pending,payload)
        job=subprocess.check_output(command,env=env,universal_newlines=True,timeout=10).strip().split(';')[0]
        if not job.isdigit():raise RuntimeError('Ambiguous source worker submission; inspect pending_submission.json')
        payload['job']=job;atomic(branch/'logs'/('finish_014_'+job+'.json'),payload)
        tmp=batch/'array_job_id.txt.online_tmp'
        with tmp.open('w') as f:f.write(job+'\n');f.flush();os.fsync(f.fileno())
        os.replace(str(tmp),str(batch/'array_job_id.txt'));pending.unlink()
        print('Finishing original 014: submitted %d %s tasks in array %s; existing outputs reused.'%(count,action,job),flush=True)
        return False

def verify_inputs(branch):
    for rel,digest in read(branch/'frozen_inputs.json').items():
        if cached_sha(branch/rel)!=digest: raise RuntimeError('Frozen branch input changed: '+rel)

def stamp(path):
    s=Path(path).stat();return dict(size=s.st_size,mtime_ns=s.st_mtime_ns)

def is_done(branch,d,gene):
    f=d/'receipts'/(gene+'.json')
    if not f.exists(): return False
    r=read(f)
    if r['gene']!=gene or r['input_sha256']!=cached_sha(d/'input.rds'): raise RuntimeError('Foreign completion receipt')
    for name,record in r['files'].items():
        p=d/name
        if not p.exists() or stamp(p)!=record['stamp']: raise RuntimeError('Completed output changed: '+str(p))
    return True

def active_assignments(branch,rows,stage):
    cfg=read(branch/'config.json');worker_name=names(cfg)[1]
    own=[r for r in rows if r['name']==worker_name]
    jobs={};active=set();attempts={}
    for p in sorted((branch/'tasks').glob('task_*/submission.json')):
        sub=read(p);task=sub['task'];jobs[sub['job']]=task
        if task['stage']==stage:
            for g in task['genes']: attempts[g]=attempts.get(g,0)+1
    for r in own:
        if r['array'] not in jobs: raise RuntimeError('Unregistered branch worker '+r['job'])
        task=jobs[r['array']]
        if task['stage']!=stage: raise RuntimeError('Worker still active for another prior/batch')
        k=int(r['job'].rsplit('_',1)[1])-1
        if k<0 or k>=len(task['genes']): raise RuntimeError('Invalid array task')
        g=task['genes'][k]
        if g in active: raise RuntimeError('Duplicate gene assignment')
        active.add(g)
    return active,attempts,len(own)

def retire(branch,current,stages):
    if legacy_jobs(queue()):raise RuntimeError('Legacy slider jobs reappeared; cleanup is blocked')
    root=result_root(branch);iteration=int(current['iteration'])
    source_iteration=int(read(branch/'config.json')['source_iteration'])
    first_online=source_iteration+1
    if iteration<first_online:raise RuntimeError('Cannot retire before the first online scan')
    if iteration>first_online:
        first_pass=[int(row['stage']) for row in stages if int(row['iteration'])==first_online]
        progress=int(csvrows(branch/'progress.csv')[0]['stage'])
        if not first_pass or progress<max(first_pass):
            raise RuntimeError('Cannot release source warm starts before the entire first online pass commits')
        for s in (branch/'seed').glob('snapshot_*'): release_results(branch,s/'results')
    for i in range(iteration-1):
        d=root/('iteration_%03d'%i)
        if not d.exists():continue
        if i>=first_online:assert_owned_iteration(branch,i)
        results=d/'results'
        if results.exists():release_results(root,results)
        marker=d/'online_retirement.json'
        if not marker.exists():atomic(marker,dict(before_iteration=iteration,history_preserved=True,workflow=str(branch)))

def control(branch):
    verify_code(branch/'code')
    if (branch/'DONE.json').exists(): print('Experiment complete; no more jobs.');return
    if (branch/'ERROR.json').exists(): raise RuntimeError('Stopped experiment: inspect ERROR.json before --resume')
    current=os.environ.get('SLURM_JOB_ID','')
    if not current.isdigit(): raise RuntimeError('Controller must run under Slurm')
    with lock(branch/'.controller.lock'):
        if (branch/'pending_submission.json').exists(): raise RuntimeError('Unresolved submission; inspect pending_submission.json')
        # A successor exists before potentially long snapshot/aggregation work.
        # TIMEOUT/OOM can therefore resume from durable transaction boundaries.
        rows=submission_queue(branch,1,20*3600)
        submit(branch,'controller',rows,dependency=current)
        if not finish_source(branch):return
        if seed_run(branch) is False: return
        if legacy_jobs(queue()):raise RuntimeError('Legacy slider jobs reappeared after the handoff; no overlapping work is allowed')
        stages=csvrows(branch/'stages.csv')
        progress=int(csvrows(branch/'progress.csv')[0]['stage'])
        if progress==59:
            atomic(branch/'DONE.json',dict(stages=59,minibatch_passes=6,full_steps=5,time=time.time()))
            print('Completed six minibatch passes and five full EM steps.',flush=True);return
        row=stages[progress];d=branch/row['directory']
        # Finalize a commit interrupted before its state pointer/history was saved.
        if (d/'commit.rds').exists():
            rrun(branch,'commit',branch,d);rrun(branch,'history',branch);return
        assert_owned_iteration(branch,int(row['iteration']))
        rrun(branch,'prepare',branch)
        rrun(branch,'history',branch)
        retire(branch,row,stages)
        genes=(branch/row['genes_file']).read_text().splitlines()
        rows=queue();active,attempts,nworkers=active_assignments(branch,rows,int(row['stage']))
        done={g for g in genes if is_done(branch,d,g)}
        print('Stage %s %s pass %s batch %s: %d/%d validated; %d active workers.'%
              (row['stage'],row['phase'],row['pass'],row['batch'],len(done),len(genes),nworkers),flush=True)
        if len(done)==len(genes) and not active:
            rrun(branch,'commit',branch,d);rrun(branch,'history',branch)
            if int(row['stage'])==59:
                atomic(branch/'DONE.json',dict(stages=59,minibatch_passes=6,full_steps=5,time=time.time()))
            return
        pending=[g for g in genes if g not in done and g not in active]
        exhausted=[g for g in pending if attempts.get(g,0)>=3]
        if exhausted: raise RuntimeError('Three worker attempts exhausted; no genes skipped: '+','.join(exhausted[:20]))
        if not pending: return
        rows=submission_queue(branch,1,20*3600)
        active,attempts,nworkers=active_assignments(branch,rows,int(row['stage']))
        pending=[g for g in genes if g not in active and not is_done(branch,d,g)]
        exhausted=[g for g in pending if attempts.get(g,0)>=3]
        if exhausted: raise RuntimeError('Three worker attempts exhausted; no genes skipped: '+','.join(exhausted[:20]))
        cfg=read(branch/'config.json')
        capacity=min(cfg['workers']-nworkers,submission_budget(rows,cfg)-len(rows))
        if capacity<=0: return
        selected=pending[:capacity]
        if not selected:return
        tasks=branch/'tasks';tasks.mkdir(exist_ok=True)
        td=tasks/('task_%06d'%(len(list(tasks.glob('task_*')))+1));td.mkdir()
        task=dict(directory=str(td),stage=int(row['stage']),stage_dir=str(d),genes=selected,
                  input_sha256=sha(d/'input.rds'))
        atomic(td/'task.json',task)
        submit(branch,'workers',rows,task=task)

def worker(branch,taskfile,index):
    verify_code(branch/'code');verify_inputs(branch);task=read(taskfile);k=int(index)-1
    if k<0 or k>=len(task['genes']): raise RuntimeError('Unknown task index')
    gene=task['genes'][k];d=Path(task['stage_dir'])
    if result_root(branch).resolve() not in d.resolve().parents or sha(d/'input.rds')!=task['input_sha256']:
        raise RuntimeError('Foreign/changed stage input')
    locks=branch/'gene_locks';locks.mkdir(exist_ok=True)
    with lock(locks/('%03d_%s.lock'%(task['stage'],gene))):
        if is_done(branch,d,gene): return
        atomic(d/'attempts'/(gene+'_'+os.environ.get('SLURM_JOB_ID','local')+'.json'),
               dict(gene=gene,stage=task['stage'],started=time.time()))
        rrun(branch,'worker',branch,d,gene)
        success=csvrows(d/'success'/(gene+'.csv'))[0]
        if success['gene']!=gene or int(success['stage'])!=task['stage']: raise RuntimeError('Wrong R receipt')
        paths=[success['result_path'],'summaries/'+gene+'.rds','success/'+gene+'.csv']
        atomic(d/'receipts'/(gene+'.json'),dict(gene=gene,input_sha256=task['input_sha256'],
            files={p:dict(stamp=stamp(d/p)) for p in paths},completed=time.time()))

def launch(a):
    bundle=Path(__file__).resolve().parent;verify_code(bundle)
    source=Path(a.source).resolve();branch=Path(a.branch).resolve()
    if branch!=source/'workflows/slider_online_015_025':raise RuntimeError('Workflow must live inside PROJECT/workflows/slider_online_015_025')
    identity=dict(schema=3,source=str(source),branch=str(branch),source_iteration=14,poll_minutes=1,rho=0.3,passes=6,batches=9,full_steps=5,
                  seed=20261104,workers=a.workers,id=hashlib.sha256(str(branch).encode()).hexdigest()[:6],
                  bundle_sha256=sha(bundle/'SHA256SUMS'))
    if branch.exists():
        if not a.resume: raise RuntimeError('Branch already exists; use --resume after inspecting its status')
        if not (branch/'config.json').exists() or read(branch/'config.json')!=identity:
            raise RuntimeError('Existing branch has different configuration')
    else:
        if a.resume: raise RuntimeError('No existing branch to resume')
        branch.mkdir(parents=True);atomic(branch/'config.json',identity)
    if not (branch/'code').exists(): shutil.copytree(str(bundle),str(branch/'code'),ignore=shutil.ignore_patterns('__pycache__'))
    verify_code(branch/'code')
    (branch/'logs').mkdir(exist_ok=True)
    with lock(branch/'.launch.lock'):
        if (branch/'ERROR.json').exists():
            raise RuntimeError('Online workflow has ERROR.json; inspect the failure before retrying')
        rows=queue();ctl,_=names(identity)
        if any(r['name']==ctl for r in rows): raise RuntimeError('A branch controller is already queued/running')
        if (branch/'DONE.json').exists(): print('Already complete.');return
        if (branch/'pending_submission.json').exists(): raise RuntimeError('Inspect ambiguous pending_submission.json before retrying')
        if not (branch/'HANDOFF.json').exists():
            if not a.handoff:raise RuntimeError('Use --handoff to transfer the pending legacy controller and finish 014')
            import handoff
            handoff.begin(branch,sys.modules[__name__])
        rows=submission_queue(branch,1,24*3600)
        submit(branch,'controller',rows)
    print('Online iterations 015..020 and full iterations 021..025 in '+str(result_root(branch)))

def status(branch):
    cfg=read(branch/'config.json')
    print('Self-contained workflow: '+str(branch))
    print('Results: '+str(result_root(branch)))
    print('Finish full 014; online 015..020; full 021..025. Refill eligibility every minute.')
    print('rho=0.3; six passes x nine batches; five full EM steps; worker cap='+str(cfg['workers']))
    ready=branch/'seed/READY.json'
    if ready.exists():
        x=read(ready);print('Source checkpoint: iteration_%03d; %d genes; %d tissue fits'%
                           (x['source_iteration'],x['genes'],x['tissue_fits']))
        print('Preserved source history: '+str(Path(x['snapshot'])/'source_prior_history.csv'))
    if (branch/'progress.csv').exists():
        progress=int(csvrows(branch/'progress.csv')[0]['stage'])
        print('Committed stages: %d/59'%progress)
        if progress<59:
            row=csvrows(branch/'stages.csv')[progress];d=branch/row['directory']
            genes=(branch/row['genes_file']).read_text().splitlines()
            completed=len(list((d/'receipts').glob('*.json')))
            print('Current stage: '+row['directory']+'; saved validation receipts: %d/%d'%(completed,len(genes)))
    for name in ('ERROR.json','pending_submission.json','DONE.json'):
        if (branch/name).exists():print(name+': '+(branch/name).read_text())
    if (branch/'last_controller.json').exists():
        job=read(branch/'last_controller.json')['job'];print('Last submitted controller: '+job)
    rows=queue();own=[r for r in rows if r['name'] in names(cfg)]
    print('Branch jobs: %d running, %d pending; all-user occupied slots: %d/300'%
          (sum(r['state']=='RUNNING' for r in own),sum(r['state']=='PENDING' for r in own),len(rows)))
    logs=sorted((branch/'logs').glob('control_*.err'),key=lambda p:p.stat().st_mtime,reverse=True)
    for p in logs[:1]:
        print('Latest controller error log: '+str(p));print('\n'.join(p.read_text(errors='replace').splitlines()[-15:]))

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='mode')
    l=sub.add_parser('launch');l.add_argument('--source',default=SOURCE);l.add_argument('--branch',default=BRANCH)
    l.add_argument('--workers',type=int,default=298);l.add_argument('--resume',action='store_true');l.add_argument('--handoff',action='store_true')
    c=sub.add_parser('control');c.add_argument('branch',type=Path)
    s=sub.add_parser('status');s.add_argument('branch',type=Path,nargs='?',default=Path(BRANCH))
    w=sub.add_parser('worker');w.add_argument('branch',type=Path);w.add_argument('taskfile');w.add_argument('index')
    a=p.parse_args()
    if a.mode=='launch':
        if not 1<=a.workers<=298: p.error('--workers must be 1..298')
        launch(a)
    elif a.mode=='worker': worker(a.branch.resolve(),a.taskfile,a.index)
    elif a.mode=='status': status(a.branch.resolve())
    elif a.mode=='control':
        branch=a.branch.resolve()
        try:control(branch)
        except Exception as e:
            if not (branch/'ERROR.json').exists():
                atomic(branch/'ERROR.json',dict(error=str(e),time=time.time(),job=os.environ.get('SLURM_JOB_ID')))
            raise
    else:p.error('Choose launch/control/worker')

if __name__=='__main__':
    try:main()
    except Exception as e: print('STOP: '+str(e),file=sys.stderr);sys.exit(1)
