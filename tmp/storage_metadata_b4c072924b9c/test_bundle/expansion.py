"""Separate expanded slider branch; reuse original completed E-step at the same priors."""
import argparse,json,os,shutil,subprocess,sys,time
from pathlib import Path
import storage

BUNDLE=Path(__file__).resolve().parent

def write(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='\n') as f: f.write(text.replace('\r\n','\n'))

def make_code(branch):
    base=BUNDLE/'baseline/script/scan_tissue_attempt'
    for source in base.glob('*.R'):
        text=source.read_text()
        if source.name=='slide_prior_recovery.R':
            old='    out <- readRDS(file)\n    sre_validate(out, ctx, previous)'
            new='''    out <- readRDS(file)
    if (ctx$iteration == 0L && !identical(as.integer(attr(out,"em_iteration")),0L)) {
      seed <- readRDS(file.path(ctx$project,"results_slide_prior_em/seed_provenance.rds"))
      expected <- seed$result_md5[gene]
      if (length(expected)!=1L || is.na(expected) || expected != row$result_md5 ||
          !identical(as.integer(attr(out,"em_iteration")),as.integer(seed$source_iteration)))
        stop("Seed result does not match its original checkpoint")
      attr(out,"em_iteration") <- 0L
    }
    sre_validate(out, ctx, previous)'''
            if text.count(old)!=1: raise RuntimeError('Unexpected recovery reader')
            text=text.replace(old,new)
        text=text.replace('slide_parallel_', 'slide_expanded_')
        write(branch/'script/scan_tissue_attempt'/source.name,text)
    for source in (BUNDLE/'script/scan_tissue_attempt').glob('*.R'):
        text=source.read_text()
        if source.name=='workhorse_slide_prior_em.R': text+='\n'+(BUNDLE/'expanded_worker_wrapper.R').read_text()
        write(branch/'script/scan_tissue_attempt'/source.name,text)
    for name in ('resume_slide_prior_em','slide_prior_parallel_worker'):
        text=(BUNDLE/'patch/job'/name).read_text().replace('slide_parallel_','slide_expanded_')
        if name=='resume_slide_prior_em':
            needle='export SUSIE_MIX_PROJECT_DIR="$project_dir"'
            text=text.replace(needle,needle+'\npython3 "$project_dir/verify_expanded.py" "$project_dir"')
        else:
            needle='project=$1; iteration=$2; batch=$3'
            text=text.replace(needle,needle+'\npython3 "$project/verify_expanded.py" "$project"')
        write(branch/'job'/name,text)
    for name in ('verify_expanded.py','renamed_gene_map.csv'):
        shutil.copyfile(BUNDLE/name,branch/name)

def initialize(a):
    source,branch=a.project.resolve(),a.branch.resolve()
    if source==branch or storage.contained(source,branch): raise RuntimeError('Branch must be a separate directory')
    storage.verify_bundle();storage.check_sources(source)
    compact=json.loads((a.state.resolve()/'COMPLETE.json').read_text())
    if compact.get('project')!=str(source) or 'saved_bytes' not in compact: raise RuntimeError('Complete storage recovery first')
    windows=json.loads((BUNDLE/'checks/genotype_summary.json').read_text())
    if windows['candidate_ids']!=680 or windows['windows_without_variants']!=0: raise RuntimeError('Genotype preflight has not passed for all candidates')
    identity=dict(source=str(source),source_iteration=a.source_iteration,bundle=str(BUNDLE))
    if (branch/'INITIALIZED.json').exists():
        if json.loads((branch/'identity.json').read_text())!=identity: raise RuntimeError('Branch identity differs')
        print('Branch already initialized; use start.',branch);return
    with storage.maintenance_lock(source):
        if branch.exists():
            if not (branch/'identity.json').exists() or json.loads((branch/'identity.json').read_text())!=identity:
                raise RuntimeError('Existing branch is unrelated; it was not overwritten')
        else: branch.mkdir(parents=True)
        storage.atomic_json(branch/'identity.json',identity)
        make_code(branch)
        subprocess.check_call([os.environ.get('RSCRIPT','Rscript'),'--vanilla',str(BUNDLE/'expanded.R'),
                               str(BUNDLE),str(source),str(branch),str(a.source_iteration)])
        dependency=source/'results_slide_prior_em/storage_dependencies.json'
        deps=json.loads(dependency.read_text()) if dependency.exists() else []
        row=dict(branch=str(branch),source_iteration=a.source_iteration)
        if row not in deps: deps.append(row)
        storage.atomic_json(dependency,deps)
        paths=list((branch/'script').rglob('*.R'))+list((branch/'job').iterdir())+[branch/'renamed_gene_map.csv',branch/'verify_expanded.py']
        paths+=[branch/'results_slide_prior_em/seed_provenance.rds']
        storage.atomic_json(branch/'frozen_expansion.json',dict(hashes={p.relative_to(branch).as_posix():storage.hashes(p)['sha256'] for p in paths}))
        storage.atomic_json(branch/'INITIALIZED.json',dict(**identity,initialized=time.time()))
        print('Expanded slider initialized. Original results and history were not changed:',branch)

def start(a):
    branch=a.branch.resolve()
    subprocess.check_call([sys.executable,str(branch/'verify_expanded.py'),str(branch)])
    with storage.maintenance_lock(branch):
        env=os.environ.copy()
        env.update(SUSIE_MIX_PROJECT_DIR=str(branch),SUSIE_SLIDE_TARGET_ITERATION=str(a.target),
                   SUSIE_SLIDE_WORKERS='298',SUSIE_SLIDE_POLL_MINUTES='5',SUSIE_SLIDE_WAIT_FOR_SIMULATIONS='0',
                   SUSIE_SLIDE_MAX_ITER='5000',SUSIE_SLIDE_RECORD_ELBO='0',
                   SUSIE_SLIDE_EXCLUDE_NODES='midway2-0003,midway2-0011')
        dependency=['--dependency=afterany:'+os.environ['SLURM_JOB_ID']] if os.environ.get('SLURM_JOB_ID') else ['--begin=now+1minutes']
        subprocess.check_call(['sbatch']+dependency+['--exclude=midway2-0003,midway2-0011','--chdir='+str(branch/'job'),
             '--output='+str(branch/'job/slide_expanded_control_%j.out'),
             '--error='+str(branch/'job/slide_expanded_control_%j.err'),str(branch/'job/resume_slide_prior_em'),'latest'],env=env)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='mode')
    q=sub.add_parser('init');q.add_argument('project',type=Path);q.add_argument('state',type=Path);q.add_argument('branch',type=Path);q.add_argument('--source-iteration',type=int,default=5)
    q=sub.add_parser('start');q.add_argument('branch',type=Path);q.add_argument('--target',type=int,default=15)
    a=p.parse_args()
    try:
        if a.mode=='init': initialize(a)
        elif a.mode=='start': start(a)
        else: p.error('Choose init or start')
    except Exception as e:
        print('STOP: '+str(e),file=sys.stderr);sys.exit(1)
