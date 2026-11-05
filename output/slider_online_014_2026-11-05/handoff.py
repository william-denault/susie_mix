"""Stop only the pending continuation; keep every fitting job running."""
import os, re, subprocess, time
from pathlib import Path

def boundary(project, mb):
    root=project/'results_slide_prior_em'
    later=[p.name for p in root.glob('iteration_*') if re.fullmatch(r'iteration_[0-9]+',p.name) and int(p.name[10:])>=14]
    if later:raise RuntimeError('Cannot take over 014: later iteration folders already exist: '+', '.join(sorted(later)))
    h=mb.csvrows(root/'prior_history.csv')
    if not h or max(int(float(r['iteration'])) for r in h)!=13:
        raise RuntimeError('The existing full EM must be at iteration 013')
    if not (root/'iteration_013/settings.rds').exists():raise RuntimeError('Iteration 013 is not prepared')

def inspect_controller(project, row, mb):
    if row['name']!='slide_inplace_control' or not row['job'].isdigit():
        raise RuntimeError('Unexpected legacy controller')
    output=subprocess.check_output(['scontrol','show','job','-o',row['job']],universal_newlines=True,timeout=15)
    fields=dict(re.findall(r'(?:^|\s)(\w+)=(.*?)(?=\s\w+=|$)',output.strip()))
    expected=str(project/'slider_inplace_runtime_2026-10-29/job/resume_slide_prior_em')
    if fields.get('JobState')!='PENDING' or fields.get('JobName')!='slide_inplace_control' or \
       fields.get('Command','').split(' ')[0]!=expected or \
       fields.get('UserId','').split('(')[0]!=os.environ['USER']:
        raise RuntimeError('Refusing to cancel a job that is not this project\'s pending controller')
    return fields

def submit_finisher(branch, mb):
    cfg=mb.read(branch/'config.json');project=Path(cfg['source']);record=branch/'finish_013_submission.json'
    if record.exists():return mb.read(record)['job']
    rows=mb.submission_queue(branch,2)
    if mb.legacy_jobs(rows) and any(r['name']=='slide_inplace_control' for r in rows):
        raise RuntimeError('A legacy controller appeared before the target-013 replacement')
    pending=branch/'pending_submission.json'
    if pending.exists():raise RuntimeError('Inspect pending_submission.json before retrying')
    # Export variables through the subprocess environment so the comma in the
    # original controller's excluded-node list cannot split --export entries.
    env=os.environ.copy();env.update(SUSIE_MIX_PROJECT_DIR=str(project),SUSIE_SLIDE_TARGET_ITERATION='13',
        SUSIE_SLIDE_WORKERS='298',SUSIE_SLIDE_POLL_MINUTES='5',SUSIE_SLIDE_WAIT_FOR_SIMULATIONS='0',
        SUSIE_SLIDE_MAX_ITER='5000',SUSIE_SLIDE_EXCLUDE_NODES='midway2-0003,midway2-0011',SBATCH_TIMELIMIT='23:00:00')
    cmd=['sbatch','--parsable','--export=ALL','--time=23:00:00','--begin=now+5minutes',
         '--job-name=slide_inplace_control','--exclude=midway2-0003,midway2-0011',
         '--chdir='+str(project/'slider_inplace_runtime_2026-10-29/job'),
         '--output='+str(branch/'logs/finish_013_%j.out'),'--error='+str(branch/'logs/finish_013_%j.err'),
         str(project/'slider_inplace_runtime_2026-10-29/job/resume_slide_prior_em'),'13']
    payload=dict(kind='finish_013',command=cmd,target=13,time=time.time())
    mb.atomic(pending,payload)
    job=subprocess.check_output(cmd,env=env,universal_newlines=True,timeout=10).strip().split(';')[0]
    if not job.isdigit():raise RuntimeError('Ambiguous target-013 submission')
    payload['job']=job;mb.atomic(record,payload);pending.unlink()
    print('Submitted existing full EM controller '+job+' with target 013; existing workers continue.',flush=True)
    return job

def begin(branch,mb):
    cfg=mb.read(branch/'config.json');project=Path(cfg['source']);root=project/'results_slide_prior_em'
    mb.verify_source(project);boundary(project,mb)
    journal=branch/'handoff_pending.json'
    if not journal.exists():
        deadline=time.monotonic()+600
        while True:
            rows=mb.queue();old=[r for r in rows if r['name']=='slide_inplace_control']
            unexpected=[r for r in mb.legacy_jobs(rows) if r['name'] not in ('slide_inplace_control','slide_inplace_fit','slide_inplace_validate')]
            if unexpected:raise RuntimeError('An older launcher/recovery chain is present; inspect it before the handoff')
            if any(r['name'].startswith('slide_mb_') for r in rows):
                raise RuntimeError('The earlier separate minibatch experiment has jobs; resolve that run before an in-place handoff')
            if len(old)>1:raise RuntimeError('More than one original controller exists; inspect the queue')
            if not old or mb.safe_window(rows,cfg):break
            if time.monotonic()>deadline:raise RuntimeError('No safe pending-controller window; retry --resume --handoff')
            print('Waiting for the current controller to finish; fitting jobs remain running.',flush=True);time.sleep(20)
        boundary(project,mb)
        fields=inspect_controller(project,old[0],mb) if old else None
        if old and (root/'parallel_continuation_job_id.txt').read_text().strip()!=old[0]['job']:
            raise RuntimeError('Queued controller differs from the recorded continuation')
        mb.atomic(journal,dict(original_job=old[0]['job'] if old else None,job_details=fields,time=time.time()))
    old_id=mb.read(journal)['original_job'];rows=mb.queue()
    old=[r for r in rows if r['job']==old_id]
    if old:
        if not mb.safe_window(rows,cfg):raise RuntimeError('Recorded controller is no longer safely pending; no cancellation performed')
        inspect_controller(project,old[0],mb);boundary(project,mb)
        subprocess.check_call(['scancel','--state=PENDING','--name=slide_inplace_control','--user='+os.environ['USER'],old_id],timeout=10)
        # No running gene job is ever passed to scancel.
        for attempt in range(6):
            if not any(r['job']==old_id for r in mb.queue()):break
            time.sleep(2)
        else:raise RuntimeError('Pending controller is still present; inspect the handoff before resuming')
    boundary(project,mb)
    job=submit_finisher(branch,mb)
    mb.atomic(branch/'HANDOFF.json',dict(source_iteration=13,first_online_iteration=14,last_online_iteration=19,
        last_full_iteration=24,cancelled_pending_controller=old_id,finisher_job=job,time=time.time()))
