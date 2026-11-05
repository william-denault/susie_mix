"""Repair missing seed links at stage zero; retain every prior and failed log."""
import csv, hashlib, importlib.util, json, os, re, shutil, subprocess, sys, time
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
PROJECT=Path('/project2/mstephens/wdenault/susie_mix')
TAG='seed_retention_2026-11-05'

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def raw_json(obj):return (json.dumps(obj,indent=2,sort_keys=True)+'\n').encode()
def write(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.patch_tmp')
    with tmp.open('wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
    os.replace(str(tmp),str(path))
def atomic(path,obj):write(path,raw_json(obj))
def csvrows(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))
def inside(path,root):
    path=Path(path)
    if path.is_symlink() or root.resolve() not in path.resolve().parents:raise RuntimeError('Unsafe patch path: '+str(path))
    return path
def move(source,dest,root):
    inside(source,root);inside(dest,root)
    if dest.exists():raise RuntimeError('Archive destination already exists: '+str(dest))
    dest.parent.mkdir(parents=True,exist_ok=True);os.rename(str(source),str(dest))
def verify_manifest(code,manifest):
    for line in manifest.decode().splitlines():
        expected,rel=line.split('  ',1)
        p=inside(code/rel,code)
        if digest(p)!=expected:raise RuntimeError('Code checksum changed: '+rel)
def module():
    spec=importlib.util.spec_from_file_location('seed_repair_runtime',HERE/'runtime/minibatch.py')
    mb=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mb;spec.loader.exec_module(mb);return mb
def idle(mb):
    rows=mb.queue()
    active=[r for r in rows if r['name'].startswith(('slide_online_','slide_inplace_','slide_parallel_','slide_mb_'))]
    if active:raise RuntimeError('Slider jobs are active; no repair performed: '+str(active[:4]))
def sources(expected,source_dir,pins):
    genes=set();missing=0
    if pins.is_symlink():raise RuntimeError('Seed results is a symlink')
    for row in expected:
        gene=row['gene']
        if not re.fullmatch(r'[A-Za-z0-9_.-]+',gene) or gene in genes or gene in ('.','..'):
            raise RuntimeError('Invalid/duplicate expected gene')
        genes.add(gene);source=Path(row['source']);target=pins/(gene+'.rds')
        if source.is_symlink() or source.resolve()!=(source_dir/(gene+'.rds')).resolve():
            raise RuntimeError('Unexpected source result: '+gene)
        info=source.stat()
        if info.st_size!=int(float(row['size'])) or abs(info.st_mtime-float(row['mtime']))>.001:
            raise RuntimeError('Source result changed: '+gene)
        if target.is_symlink():raise RuntimeError('Symlink in seed results: '+gene)
        if target.exists():
            if not os.path.samefile(str(source),str(target)):raise RuntimeError('Seed is not the original hard link: '+gene)
        else:missing+=1
    names={g+'.rds' for g in genes}
    if pins.exists() and any(p.name not in names for p in pins.iterdir()):
        raise RuntimeError('Unexpected seed result files')
    return missing
def restore(expected,pins):
    pins.mkdir(parents=True,exist_ok=True);created=0
    for i,row in enumerate(expected,1):
        target=pins/(row['gene']+'.rds')
        if not target.exists():os.link(row['source'],str(target),follow_symlinks=False);created+=1
        if i%1000==0:print('Restored/verified %d/%d seed links (no RDS copies).'%(i,len(expected)),flush=True)
    return created
def patch_code(branch,state,old_manifest,new_manifest):
    code=branch/'code';saved=state/'code.before'
    if code.exists():
        current=(code/'SHA256SUMS').read_bytes()
        if current==new_manifest:verify_manifest(code,new_manifest);return
        if current!=old_manifest:raise RuntimeError('Unexpected installed code manifest')
        verify_manifest(code,old_manifest)
    elif not saved.exists():raise RuntimeError('Neither installed code nor its backup exists')
    stage=state/'code.next'
    if not stage.exists():shutil.copytree(str(HERE/'runtime'),str(stage),ignore=shutil.ignore_patterns('__pycache__'))
    verify_manifest(stage,new_manifest)
    if code.exists():move(code,saved,branch)
    move(stage,code,branch)
def resume(branch,state,mb):
    mb.verify_code(branch/'code');mb.verify_inputs(branch)
    if (branch/'pending_submission.json').exists():raise RuntimeError('Ambiguous submission; inspect pending_submission.json')
    rows=mb.queue();cfg=read(branch/'config.json')
    if any(r['name']==mb.names(cfg)[0] for r in rows):
        print('A workflow controller is already queued/running; no second submission.',flush=True);return
    # The patched launcher supplies the same configuration, using its updated manifest.
    subprocess.check_call(['bash',str(branch/'code/start.sh'),'--resume','--handoff'])

def repair(project=PROJECT,expected_genes=18263):
    project=project.resolve();branch=project/'workflows/slider_online_015_025'
    inside(branch,project)
    old_manifest=(HERE/'baseline_SHA256SUMS').read_bytes();new_manifest=(HERE/'runtime/SHA256SUMS').read_bytes()
    verify_manifest(HERE,(HERE/'SHA256SUMS').read_bytes());mb=module()
    state=branch/'repairs'/TAG;state.mkdir(parents=True,exist_ok=True)
    with mb.lock(branch/'.seed_repair.lock'):
      if (state/'COMPLETE.json').exists():
        print('Repair already completed; checking continuation.',flush=True)
      else:
       with mb.lock(branch/'.launch.lock'),mb.lock(branch/'.controller.lock'):
        idle(mb)
        if (branch/'pending_submission.json').exists():raise RuntimeError('Unresolved submission; no automatic repair')
        baseline=state/'PLAN.json'
        if not baseline.exists():
            verify_manifest(branch/'code',old_manifest);mb.verify_inputs(branch)
            cfg=read(branch/'config.json')
            if cfg.get('schema')!=3 or cfg.get('source_iteration')!=14 or Path(cfg['source']).resolve()!=project or \
               cfg.get('rho')!=.3 or cfg.get('passes')!=6 or cfg.get('batches')!=9 or cfg.get('full_steps')!=5 or \
               cfg.get('bundle_sha256')!=hashlib.sha256(old_manifest).hexdigest():
                raise RuntimeError('Not the expected 015 workflow/configuration')
            error=read(branch/'ERROR.json')
            if not error['error'].startswith('Three worker attempts exhausted; no genes skipped:'):
                raise RuntimeError('Different failure; refusing to clear ERROR.json')
            if int(csvrows(branch/'progress.csv')[0]['stage'])!=0:raise RuntimeError('A prior update has already committed')
            ready=read(branch/'seed/READY.json');snapshot=Path(ready['snapshot'])
            if snapshot.resolve()!=(branch/'seed/snapshot_001').resolve() or ready['source_iteration']!=14 or ready['genes']!=expected_genes:
                raise RuntimeError('Unexpected source snapshot')
            first=project/'results_slide_prior_em/iteration_015/batches/batch_01'
            results=first.parents[1]/'results'
            if any(results.iterdir()) or any((first/'success').iterdir()) or \
               ((first/'receipts').exists() and any((first/'receipts').iterdir())) or (first/'commit.rds').exists():
                raise RuntimeError('015 already contains results/receipts; use a broader recovery')
            for task in (branch/'tasks').glob('task_*'):
                if task.is_symlink() or int(read(task/'task.json')['stage'])!=1:raise RuntimeError('Unexpected task history')
            expected=csvrows(snapshot/'expected.csv')
            if len(expected)!=expected_genes:raise RuntimeError('Wrong seed cohort')
            missing=sources(expected,project/'results_slide_prior_em/iteration_014/results',snapshot/'results')
            if not missing:raise RuntimeError('No seed links are missing; this failure needs a different diagnosis')
            subprocess.check_call([os.environ.get('RSCRIPT','Rscript'),'--vanilla',str(HERE/'check_state.R'),str(branch),str(expected_genes)])
            before=state/'before';before.mkdir(exist_ok=True)
            for rel in ('config.json','frozen_inputs.json','ERROR.json'):
                shutil.copyfile(str(branch/rel),str(before/rel))
            keep=[branch/'state.rds',branch/'model.rds',branch/'stages.csv',branch/'initial_priors.csv',
                  branch/'original_prior_history.csv',branch/'progress.csv',first/'input.rds',
                  project/'results_slide_prior_em/prior_history.csv',
                  project/'results_slide_prior_em/prior_history_long.csv']
            cfg['bundle_sha256']=hashlib.sha256(new_manifest).hexdigest()
            new_config=raw_json(cfg);new_frozen=read(branch/'frozen_inputs.json')
            new_frozen['config.json']=hashlib.sha256(new_config).hexdigest()
            write(state/'config.after.json',new_config);atomic(state/'frozen.after.json',new_frozen)
            atomic(baseline,dict(project=str(project),branch=str(branch),snapshot=str(snapshot),missing_before=missing,
                old_manifest_sha256=hashlib.sha256(old_manifest).hexdigest(),
                new_manifest_sha256=hashlib.sha256(new_manifest).hexdigest(),
                unchanged={str(p):digest(p) for p in keep},time=time.time()))
        plan=read(baseline)
        if plan['project']!=str(project) or plan['branch']!=str(branch) or \
           plan['new_manifest_sha256']!=hashlib.sha256(new_manifest).hexdigest():raise RuntimeError('Foreign repair transaction')
        for name,sha in plan['unchanged'].items():
            if digest(name)!=sha:raise RuntimeError('State/history changed during repair: '+name)
        snapshot=Path(plan['snapshot']);expected=csvrows(snapshot/'expected.csv')
        sources(expected,project/'results_slide_prior_em/iteration_014/results',snapshot/'results')
        restore(expected,snapshot/'results')
        if sources(expected,project/'results_slide_prior_em/iteration_014/results',snapshot/'results'):
            raise RuntimeError('Some seed links are still missing')
        idle(mb);patch_code(branch,state,old_manifest,new_manifest)
        for rel,after in [('config.json','config.after.json'),('frozen_inputs.json','frozen.after.json')]:
            now=(branch/rel).read_bytes();old=(state/'before'/rel).read_bytes();new=(state/after).read_bytes()
            if now not in (old,new):raise RuntimeError('Changed configuration during repair: '+rel)
            write(branch/rel,new)
        mb.HASH_CACHE.clear();mb.verify_code(branch/'code');mb.verify_inputs(branch)
        tasks=branch/'tasks';archived=state/'failed_tasks'
        if not archived.exists():move(tasks,archived,branch)
        elif tasks.exists() and any(tasks.iterdir()):raise RuntimeError('New tasks appeared during repair')
        tasks.mkdir(exist_ok=True)
        error=branch/'ERROR.json'
        if error.exists():
            if error.read_bytes()!=(state/'before/ERROR.json').read_bytes():raise RuntimeError('A different error appeared during repair')
            move(error,state/'ERROR.before_repair.json',branch)
        for name,sha in plan['unchanged'].items():
            if digest(name)!=sha:raise RuntimeError('Repair changed frozen analysis state: '+name)
        atomic(state/'COMPLETE.json',dict(seed_links=len(expected),original_iteration=14,resume_iteration=15,
            committed_updates=0,priors_and_counts_unchanged=True,failed_tasks_archived=True,time=time.time()))
        print('Repair verified: %d seed links restored; priors/counts unchanged; failed attempts archived.'%len(expected),flush=True)
      resume(branch,state,mb)

if __name__=='__main__':
    try:repair()
    except Exception as exc:print('STOP: '+str(exc),file=sys.stderr);sys.exit(1)
