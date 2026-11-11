"""Verified installation helpers; no fitting or scheduling at import time."""
import hashlib, importlib.util, json, os, re, shutil, subprocess, sys, time
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent

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
    files.update(branch.glob('publication_*.rds'))
    results={};validated=[]
    stages=mb.csvrows(branch/'stages.csv')
    if [(int(r['stage']),int(r['iteration']),int(r['batch'])) for r in stages[:6]]!=[(i,15,i) for i in range(1,7)]:
        raise RuntimeError('Unexpected batch plan')
    # Five immutable commits exist; all six batches are fully validated.
    for row in stages[6:]:
        if (branch/row['directory']/'input.rds').exists():raise RuntimeError('A later batch is already prepared')
    for index,row in enumerate(stages[:6],1):
        d=branch/row['directory'];gene_list=(branch/row['genes_file']).read_text().splitlines()
        if d.resolve()!=(root/('iteration_015/batches/batch_%02d'%index)).resolve():
            raise RuntimeError('Unexpected stage directory')
        if (d/'commit.rds').exists()!=(index<=5):raise RuntimeError('Unexpected committed batch boundary')
        files.update(d.glob('*.csv'));files.add(d/'input.rds')
        files.update(d.glob('*.rds'))
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
        if done!=set(gene_list):raise RuntimeError('A saved batch has missing validated receipts')
        validated.append(len(done))
        print('Checked batch %d: %d/%d validated result receipts.'%(index,len(done),len(gene_list)),flush=True)
    actual={str(p.resolve()) for p in (root/'iteration_015/results').glob('*.rds')}
    if actual!=set(results):raise RuntimeError('015 contains unvalidated or missing result files')
    hashes={}
    for number,p in enumerate(sorted(files),1):
        hashes[str(p)]=digest(p)
        if number%5000==0:print('Protected metadata checksums: %d/%d.'%(number,len(files)),flush=True)
    return dict(hashes=hashes,outputs=results,validated_per_batch=validated)


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


