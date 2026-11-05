"""Reconcile the accepted task_000038 array without replaying EM updates."""
import hashlib, importlib.util, json, os, re, shutil, subprocess, sys, time
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
PROJECT = Path('/project2/mstephens/wdenault/susie_mix')
TAG = 'submission_ack_2026-11-05'
ARRAY = '49853903'

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def read(path): return json.loads(Path(path).read_text())
def encoded(obj): return (json.dumps(obj,indent=2,sort_keys=True)+'\n').encode()
def write(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.repair_tmp')
    with temp.open('wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
    os.replace(str(temp),str(path))
def atomic(path,obj):write(path,encoded(obj))
def inside(path,root):
    path=Path(path)
    if path.is_symlink() or root.resolve() not in path.resolve().parents:
        raise RuntimeError('Unsafe repair path: '+str(path))
    return path
def move(source,destination,root):
    inside(source,root);inside(destination,root)
    if destination.exists():raise RuntimeError('Archive already exists: '+str(destination))
    destination.parent.mkdir(parents=True,exist_ok=True)
    os.rename(str(source),str(destination))
def verify_manifest(code,manifest):
    for line in manifest.decode().splitlines():
        expected,name=line.split('  ',1)
        if digest(inside(code/name,code))!=expected:raise RuntimeError('Code changed: '+name)
def runtime():
    spec=importlib.util.spec_from_file_location('submission_repair_runtime',HERE/'runtime/minibatch.py')
    mb=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mb;spec.loader.exec_module(mb);return mb
def idle(mb):
    rows=mb.queue()
    active=[r for r in rows if r['name'].startswith(('slide_online_','slide_inplace_','slide_parallel_','slide_mb_'))
            or r['name'] in ('slide_next_gate','slide_next_setup','slide_seed_repair')]
    if active:raise RuntimeError('Slider jobs are active; leaving them unchanged: '+str(active[:4]))

def accounting(branch,mb,count):
    cfg=read(branch/'config.json');user=os.environ['USER']
    out=subprocess.check_output(['sacct','-n','-P','-X','--array','-j',ARRAY,
        '--format=JobID%40,JobName%100,User%128,State%40,ExitCode,WorkDir%1024'],
        universal_newlines=True,timeout=120)
    seen=set()
    for line in out.splitlines():
        if not line.strip():continue
        parts=[p.strip() for p in line.split('|')]
        if len(parts)!=6:raise RuntimeError('Unexpected accounting response')
        job,name,owner,status,exitcode,workdir=parts
        match=re.fullmatch(ARRAY+r'_(\d+)',job)
        if not match or int(match[1]) in seen:raise RuntimeError('Unexpected/duplicate accounting job: '+job)
        seen.add(int(match[1]))
        if name!=mb.names(cfg)[1] or owner!=user or status!='COMPLETED' or exitcode!='0:0' or \
           Path(workdir).resolve()!=branch.resolve():
            raise RuntimeError('Array identity/state differs from the reported successful run: '+line)
    if seen!=set(range(1,count+1)):
        raise RuntimeError('Accounting does not yet confirm all %d completed workers'%count)
    return out

def validate_pending(branch,mb,count):
    cfg=read(branch/'config.json');pending=read(branch/'pending_submission.json');task=pending.get('task',{})
    td=branch/'tasks/task_000038';stage=Path(cfg['source'])/'results_slide_prior_em/iteration_015/batches/batch_03'
    error=read(branch/'ERROR.json')
    if pending.get('kind')!='workers' or 'job' in pending or Path(task.get('directory','')).resolve()!=td.resolve() or \
       Path(task.get('stage_dir','')).resolve()!=stage.resolve() or task.get('stage')!=3:
        raise RuntimeError('Not the reported stage-3 worker submission')
    genes=task.get('genes',[])
    if len(genes)!=count or len(set(genes))!=count or task!=read(td/'task.json'):
        raise RuntimeError('Pending task content changed')
    command=['sbatch','--parsable','--export=ALL','--time=23:00:00',
        '--exclude=midway2-0003,midway2-0011','--chdir='+str(branch),
        '--job-name='+mb.names(cfg)[1],'--output='+str(td/'worker_%A_%a.out'),
        '--error='+str(td/'worker_%A_%a.err'),'--array=1-'+str(count),
        str(branch/'code/worker_job'),str(branch),str(td/'task.json')]
    if pending.get('command')!=command or error.get('job')!='49853687' or \
       error.get('error')!=str(subprocess.TimeoutExpired(command,10)):
        raise RuntimeError('Different failure/command; refusing to clear it')
    if mb.sha(stage/'input.rds')!=task['input_sha256']:raise RuntimeError('Batch-3 input changed')
    if (stage/'commit.rds').exists():raise RuntimeError('Batch 3 is already committed')
    observed=set()
    for p in td.glob('worker_*.err'):
        match=re.fullmatch(r'worker_(\d+)_(\d+)\.err',p.name)
        if not match or match[1]!=ARRAY:raise RuntimeError('Unexpected job logs in task_000038')
        observed.add(int(match[2]))
    if observed!=set(range(1,count+1)):raise RuntimeError('Missing worker logs for accepted array '+ARRAY)
    for gene in genes:
        if not mb.is_done(branch,stage,gene):raise RuntimeError('Accepted worker lacks a validated receipt: '+gene)
    return pending

def protected_state(branch,mb):
    """Hash metadata; only stat large RDS results, using existing receipts."""
    root=Path(read(branch/'config.json')['source'])/'results_slide_prior_em'
    files={branch/n for n in ('state.rds','model.rds','progress.csv','stages.csv','prior_history.csv',
        'initial_priors.csv','original_prior_history.csv','latest_priors.csv','HANDOFF.json','seed/READY.json')}
    files.update(root/n for n in ('prior_history.csv','prior_history_long.csv','online_latest_priors.csv'))
    files.update((root/'iteration_015').glob('*.csv'))
    files.add(root/'iteration_015/ONLINE_OWNER.json')
    files.update((branch/'tasks').glob('task_*/task.json'))
    files.update((branch/'tasks').glob('task_*/submission.json'))
    results={};validated=[]
    stages=mb.csvrows(branch/'stages.csv')
    if [(int(r['stage']),int(r['iteration']),int(r['batch'])) for r in stages[:3]]!=[(1,15,1),(2,15,2),(3,15,3)]:
        raise RuntimeError('Unexpected batch plan')
    # The first two immutable commits must already exist. Batch 3 must not.
    for index,row in enumerate(stages[:3],1):
        d=branch/row['directory'];gene_list=(branch/row['genes_file']).read_text().splitlines()
        if d.resolve()!=(root/('iteration_015/batches/batch_%02d'%index)).resolve():
            raise RuntimeError('Unexpected stage directory')
        if (d/'commit.rds').exists()!=(index<=2):raise RuntimeError('Unexpected committed batch boundary')
        files.update(d.glob('*.csv'));files.add(d/'input.rds')
        if index<=2:files.add(d/'commit.rds')
        receipts=list((d/'receipts').glob('*.json'))
        done=set()
        for p in receipts:
            gene=p.stem
            if gene not in gene_list or not mb.is_done(branch,d,gene):raise RuntimeError('Invalid receipt: '+str(p))
            done.add(gene);files.add(p)
            for name in read(p)['files']:
                target=(d/name).resolve();inside(target,root/'iteration_015')
                if target.parent==(root/'iteration_015/results').resolve():results[str(target)]=mb.stamp(target)
                else:files.add(target)
        if index<=2 and done!=set(gene_list):raise RuntimeError('A committed batch has missing receipts')
        validated.append(len(done))
    actual={str(p.resolve()) for p in (root/'iteration_015/results').glob('*.rds')}
    if actual!=set(results):raise RuntimeError('015 contains unvalidated or missing result files')
    files.discard(branch/'tasks/task_000038/submission.json')
    return dict(hashes={str(p):digest(p) for p in sorted(files)},outputs=results,
                validated_per_batch=validated)

def check_protected(plan,mb):
    for p,expected in plan['protected']['hashes'].items():
        if digest(p)!=expected:raise RuntimeError('Analysis metadata changed: '+p)
    for p,expected in plan['protected']['outputs'].items():
        if mb.stamp(p)!=expected:raise RuntimeError('A saved fit changed: '+p)

def patch_code(branch,state,old,new):
    code=branch/'code';backup=state/'code.before'
    if code.exists():
        current=(code/'SHA256SUMS').read_bytes()
        if current==new:verify_manifest(code,new);return
        if current!=old:raise RuntimeError('Unexpected installed code')
        verify_manifest(code,old)
    elif not backup.exists():raise RuntimeError('Installed code and backup both missing')
    stage=state/'code.next'
    if not stage.exists():shutil.copytree(str(HERE/'runtime'),str(stage),ignore=shutil.ignore_patterns('__pycache__'))
    verify_manifest(stage,new)
    if code.exists():move(code,backup,branch)
    move(stage,code,branch)

def resume(branch,mb):
    mb.verify_code(branch/'code');mb.verify_inputs(branch)
    if (branch/'ERROR.json').exists() or (branch/'pending_submission.json').exists():
        raise RuntimeError('A new unresolved failure is present; not resuming')
    if any(r['name']==mb.names(read(branch/'config.json'))[0] for r in mb.queue()):
        print('Controller already queued/running; no duplicate submission.',flush=True);return
    subprocess.check_call(['bash',str(branch/'code/start.sh'),'--resume','--handoff'])

def repair(project=PROJECT,expected_workers=246):
    project=project.resolve();branch=project/'workflows/slider_online_015_025';inside(branch,project)
    mb=runtime();verify_manifest(HERE,(HERE/'SHA256SUMS').read_bytes())
    old=(HERE/'baseline_SHA256SUMS').read_bytes();new=(HERE/'runtime/SHA256SUMS').read_bytes()
    state=branch/'repairs'/TAG;state.mkdir(parents=True,exist_ok=True)
    with mb.lock(branch/'.submission_repair.lock'):
      if (state/'COMPLETE.json').exists():
        print('Submission repair already complete; checking continuation.',flush=True)
      else:
       with mb.lock(branch/'.launch.lock'),mb.lock(branch/'.controller.lock'):
        idle(mb);planfile=state/'PLAN.json'
        if not planfile.exists():
            print('Checking code, saved updates, completed workers and receipts.',flush=True)
            verify_manifest(branch/'code',old);mb.verify_inputs(branch)
            cfg=read(branch/'config.json')
            expected={'schema':3,'source_iteration':14,'rho':.3,'passes':6,'batches':9,'full_steps':5,'workers':298,'poll_minutes':1}
            if any(cfg.get(k)!=v for k,v in expected.items()) or Path(cfg['source']).resolve()!=project or \
               cfg.get('bundle_sha256')!=hashlib.sha256(old).hexdigest():raise RuntimeError('Unexpected installed workflow')
            if int(mb.csvrows(branch/'progress.csv')[0]['stage'])!=2:raise RuntimeError('Expected exactly two committed batches')
            pending=validate_pending(branch,mb,expected_workers)
            evidence=accounting(branch,mb,expected_workers)
            protected=protected_state(branch,mb)
            before=state/'before';before.mkdir(exist_ok=True)
            for rel in ('config.json','frozen_inputs.json','ERROR.json','pending_submission.json'):
                shutil.copyfile(str(branch/rel),str(before/rel))
            write(state/'accounting.txt',evidence.encode())
            cfg['bundle_sha256']=hashlib.sha256(new).hexdigest();cfgbytes=encoded(cfg)
            frozen=read(branch/'frozen_inputs.json');frozen['config.json']=hashlib.sha256(cfgbytes).hexdigest()
            write(state/'config.after.json',cfgbytes);atomic(state/'frozen.after.json',frozen)
            pending['job']=ARRAY
            pending['recovered_acknowledgement']=dict(time=time.time(),evidence='accounting.txt and worker logs',repair=TAG)
            atomic(state/'submission.after.json',pending)
            atomic(planfile,dict(project=str(project),branch=str(branch),array=ARRAY,workers=expected_workers,
                new_manifest=hashlib.sha256(new).hexdigest(),protected=protected,time=time.time()))
        plan=read(planfile)
        if plan['project']!=str(project) or plan['branch']!=str(branch) or plan['array']!=ARRAY or \
           plan['new_manifest']!=hashlib.sha256(new).hexdigest():raise RuntimeError('Foreign repair transaction')
        check_protected(plan,mb);idle(mb)
        # Reconfirm termination even when recovering an interrupted installation.
        accounting(branch,mb,plan['workers'])
        patch_code(branch,state,old,new)
        for rel,after in [('config.json','config.after.json'),('frozen_inputs.json','frozen.after.json')]:
            original=(state/'before'/rel).read_bytes();updated=(state/after).read_bytes()
            if (branch/rel).read_bytes() not in (original,updated):raise RuntimeError('Changed configuration during repair')
            write(branch/rel,updated)
        receipt=branch/'tasks/task_000038/submission.json';data=(state/'submission.after.json').read_bytes()
        if receipt.exists() and receipt.read_bytes()!=data:raise RuntimeError('Different submission receipt already exists')
        write(receipt,data)
        for rel in ('pending_submission.json','ERROR.json'):
            current=branch/rel;archive=state/(rel+'.before_repair')
            expected=(state/'before'/rel).read_bytes()
            if current.exists():
                if current.read_bytes()!=expected:raise RuntimeError('Different failure appeared during repair')
                move(current,archive,branch)
            elif not archive.exists() or archive.read_bytes()!=expected:raise RuntimeError('Missing transaction journal: '+rel)
        mb.HASH_CACHE.clear();mb.verify_code(branch/'code');mb.verify_inputs(branch);check_protected(plan,mb)
        n=sum(plan['protected']['validated_per_batch'])
        atomic(state/'COMPLETE.json',dict(array=ARRAY,workers=plan['workers'],preserved_validated_results=n,
            committed_stages=2,resume_iteration=15,resume_batch=3,submit_timeout_seconds=120,time=time.time()))
        print('Recovered array %s: %d completed workers; %d validated results preserved; two updates unchanged.'%
              (ARRAY,plan['workers'],n),flush=True)
        print('Resuming iteration 015, batch 3; no completed gene is resubmitted.',flush=True)
      resume(branch,mb)

if __name__=='__main__':
    try:repair()
    except Exception as exc:print('STOP: '+str(exc),file=sys.stderr);sys.exit(1)
