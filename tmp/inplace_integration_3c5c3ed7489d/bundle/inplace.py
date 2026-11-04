"""Private code, same results folder; original 005 -> expanded 006..040."""
import argparse,json,os,shutil,subprocess,sys,time
from pathlib import Path
import storage
from common import RUNTIME,META,exclusive,recovery_ready,check_dependencies
BUNDLE=Path(__file__).resolve().parent

def write(p,text):p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(text.replace('\r\n','\n').encode())
def once(text,old,new):
    if text.count(old)!=1:raise RuntimeError('Unexpected baseline: '+old)
    return text.replace(old,new)

def make_code(project):
    runtime=project/RUNTIME
    def private(text):return text.replace('script/scan_tissue_attempt',RUNTIME+'/script/scan_tissue_attempt')
    for source in (BUNDLE/'baseline/script/scan_tissue_attempt').glob('*.R'):
        text=source.read_text()
        if source.name=='slide_prior_recovery.R':
            text=once(text,'  if (!file.exists(file)) return(list(md5 = NULL, excluded = character(), chromosomes = NULL))',
              '''  expanded_file <- file.path(project,"results_slide_prior_em/inplace_expansion_2026-10-29/chromosome_scope.rds")
  if (file.exists(expanded_file) && identical(readRDS(expanded_file)$genes$gene,manifest$gene)) file <- expanded_file
  if (!file.exists(file)) return(list(md5 = NULL, excluded = character(), chromosomes = NULL))''')
            text=once(text,'sre_previous <- function(ctx, gene) {\n  if (ctx$iteration == 0L)',
              'sre_previous <- function(ctx, gene) {\n  if (si_new_at_six(ctx,gene)) return(NULL)\n  if (ctx$iteration == 0L)')
        if source.name=='slide_prior_parallel.R':
            text=once(text,'spp_previous_cohort <- function(ctx, gene) {\n  if (ctx$iteration == 0L)',
              'spp_previous_cohort <- function(ctx, gene) {\n  if (si_new_at_six(ctx,gene)) return(NULL)\n  if (ctx$iteration == 0L)')
            text=once(text,'    if (!identical(cohort(prev$audit), cohort(pooled$audit))) stop("Frozen gene/tissue cohort changed.")',
              '''    comparison <- pooled$audit
    if (ctx$iteration == 6L) {
      original_genes <- si_info(ctx$project)$original_genes
      comparison <- comparison[comparison$gene %in% original_genes,,drop=FALSE]
    }
    if (!identical(cohort(prev$audit), cohort(comparison))) stop("Frozen gene/tissue cohort changed.")''')
        if source.name=='run_slide_prior_parallel.R':
            text=once(text,'"prepare_slide_prior_em_iteration.R", "slide_prior_parallel.R"))',
                      '"prepare_slide_prior_em_iteration.R", "slide_prior_parallel.R", "cohort.R"))')
            text=once(text,'  spp_init(ctx)\n  plan <- spp_plan',
              '  spp_init(ctx)\n  source(file.path(project,"script/scan_tissue_attempt/rolling_hook.R"))\n  sr_barrier(ctx,args[5])\n  plan <- spp_plan')
            text=once(text,'    spp_init(ctx)\n    # No old worker',
              '    spp_init(ctx)\n    sr_barrier(ctx,args[5])\n    # No old worker')
        write(runtime/'script/scan_tissue_attempt'/source.name,private(text.replace('slide_parallel_','slide_inplace_')))
    for source in (BUNDLE/'script/scan_tissue_attempt').glob('*.R'):
        text=source.read_text()
        if source.name=='workhorse_slide_prior_em.R':text+='\n'+(BUNDLE/'worker_wrapper.R').read_text()
        write(runtime/'script/scan_tissue_attempt'/source.name,private(text))
    for name in ('cohort.R','rolling_hook.R'):
        write(runtime/'script/scan_tissue_attempt'/name,(BUNDLE/name).read_text())
    for name in ('resume_slide_prior_em','slide_prior_parallel_worker'):
        text=(BUNDLE/'patch/job'/name).read_text().replace('slide_parallel_','slide_inplace_')
        if name=='resume_slide_prior_em':
            shared='''if [[ -f "$project_dir/results_slide_prior_em/shared_scheduler/enabled.json" ]]; then
    module load R/4.2.0
    exec python3 "$project_dir/script/scan_tissue_attempt/slide_shared_scheduler.py" control "$project_dir"
fi'''
            text=once(text,shared,'# This private controller counts all user jobs; the old simulation handover stays unchanged.')
            needle='export SUSIE_MIX_PROJECT_DIR="$project_dir"'
            text=once(text,needle,needle+'\npython3 "$project_dir/'+RUNTIME+'/verify_inplace.py" "$project_dir"')
            text=once(text,'slide_inplace_control|slide_recover_*|',
              'slide_inplace_control|slide_parallel_*|slide_expanded_*|slide_recover_*|')
            text=text.replace('target=${SUSIE_SLIDE_TARGET_ITERATION:-20}','target=${SUSIE_SLIDE_TARGET_ITERATION:-40}')
        else:
            needle='project=$1; iteration=$2; batch=$3'
            text=once(text,needle,needle+'\npython3 "$project/'+RUNTIME+'/verify_inplace.py" "$project"')
        text=private(text).replace('$project_dir/job','$project_dir/'+RUNTIME+'/job')
        write(runtime/'job'/name,text)
    for name in ('verify_inplace.py','common.py','renamed_gene_map.csv','storage.py','maintenance.R','retention.py','retention.R'):
        shutil.copyfile(BUNDLE/name,runtime/name)

def freeze(project):
    runtime=project/RUNTIME;meta=project/META
    files=list((runtime/'script').rglob('*.R'))+list((runtime/'job').iterdir())
    files += [runtime/n for n in ('verify_inplace.py','common.py','renamed_gene_map.csv','storage.py','maintenance.R','retention.py','retention.R')]
    files += [meta/n for n in ('identity.json','provenance.rds','manifest.csv','chromosome_scope.rds','prior_history.before_expansion.csv')]
    storage.atomic_json(runtime/'frozen.json',dict(project=str(project),hashes={p.relative_to(project).as_posix():storage.hashes(p)['sha256'] for p in files}))

def initialize(project,state):
    project,state=project.resolve(),state.resolve();meta=project/META
    storage.verify_bundle();recovery_ready(project,state)
    if (meta/'INITIALIZED.json').exists():
        subprocess.check_call([sys.executable,str(project/RUNTIME/'verify_inplace.py'),str(project)]);return
    storage.check_sources(project)
    genotype=json.loads((BUNDLE/'checks/genotype_summary.json').read_text())
    if genotype['candidate_ids']!=680 or genotype['windows_without_variants']!=0:raise RuntimeError('Candidate genotype preflight is incomplete')
    identity=dict(project=str(project),recovery_state=str(state),bundle=str(BUNDLE),first=6,target=40)
    with exclusive(project):
        check_dependencies(project)
        if meta.exists():
            if not (meta/'identity.json').exists() or json.loads((meta/'identity.json').read_text())!=identity:raise RuntimeError('Existing expansion state is unrelated')
        else:
            if (project/'results_slide_prior_em/iteration_006').exists():raise RuntimeError('An unowned iteration 006 already exists')
            if (project/RUNTIME).exists():raise RuntimeError('An unowned runtime directory already exists')
            meta.mkdir(parents=True)
        storage.atomic_json(meta/'identity.json',identity)
        make_code(project)
        subprocess.check_call([os.environ.get('RSCRIPT','Rscript'),'--vanilla',str(BUNDLE/'initialize.R'),str(BUNDLE),str(project)])
        freeze(project)
        storage.atomic_json(meta/'INITIALIZED.json',dict(**identity,initialized=time.time()))
        print('Prepared iteration_006 in the EXISTING results_slide_prior_em folder; 005 fits were not rewritten.')

def start(project):
    project=project.resolve();meta=project/META
    subprocess.check_call([sys.executable,str(project/RUNTIME/'verify_inplace.py'),str(project)])
    with exclusive(project):
        pending=meta/'controller_submission.pending.json'
        if pending.exists():raise RuntimeError('Unresolved submission: '+str(pending))
        env=os.environ.copy();env.update(SUSIE_MIX_PROJECT_DIR=str(project),SUSIE_SLIDE_TARGET_ITERATION='40',
          SUSIE_SLIDE_WORKERS='298',SUSIE_SLIDE_POLL_MINUTES='5',SUSIE_SLIDE_WAIT_FOR_SIMULATIONS='0',
          SUSIE_SLIDE_MAX_ITER='5000',SUSIE_SLIDE_RECORD_ELBO='0',SUSIE_SLIDE_EXCLUDE_NODES='midway2-0003,midway2-0011')
        dependency=['--dependency=afterany:'+os.environ['SLURM_JOB_ID']] if os.environ.get('SLURM_JOB_ID') else ['--begin=now+1minutes']
        jobdir=project/RUNTIME/'job'
        cmd=['sbatch','--parsable']+dependency+['--exclude=midway2-0003,midway2-0011','--chdir='+str(jobdir),
          '--output='+str(jobdir/'slide_inplace_control_%j.out'),'--error='+str(jobdir/'slide_inplace_control_%j.err'),
          str(jobdir/'resume_slide_prior_em'),'latest']
        storage.atomic_json(pending,dict(command=cmd,target=40))
        job=subprocess.check_output(cmd,env=env,universal_newlines=True).strip().split(';')[0]
        if not job.isdigit():raise RuntimeError('Ambiguous controller submission')
        storage.atomic_json(meta/'last_controller.json',dict(job=job,target=40));pending.unlink();storage.sync_dir(meta)
        print('Submitted controller '+job+'; final target iteration_040, same results folder.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['init','start'])
    p.add_argument('project',type=Path);p.add_argument('state',type=Path,nargs='?');a=p.parse_args()
    try:
        if a.mode=='init':initialize(a.project,a.state)
        else:start(a.project)
    except Exception as e:print('STOP: '+str(e),file=sys.stderr);sys.exit(1)
