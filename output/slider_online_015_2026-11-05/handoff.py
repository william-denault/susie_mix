"""Hold and verify a pending legacy controller; preserve every gene worker."""
import os, re, subprocess, time
from pathlib import Path

def boundary(project, mb):
    root=project/'results_slide_prior_em'
    later=[p.name for p in root.glob('iteration_*')
           if re.fullmatch(r'iteration_[0-9]+',p.name) and int(p.name[10:])>=15]
    if later:raise RuntimeError('Cannot take over 015: later iteration folders already exist: '+', '.join(sorted(later)))
    history=mb.csvrows(root/'prior_history.csv')
    if not history or max(int(float(r['iteration'])) for r in history)!=14:
        raise RuntimeError('The existing full EM must be at iteration 014')
    if not (root/'iteration_014/settings.rds').exists():raise RuntimeError('Iteration 014 is not prepared')
    receipt=root/'rolling_cleanup_inplace/before_014/COMPLETE.json'
    if not receipt.exists() or mb.read(receipt).get('iteration')!=14:
        raise RuntimeError('Original cleanup before 014 is not certified; leave its controller in place')
    prior=project/'workflows/slider_online_014_024'
    if any((prior/p).exists() for p in ('HANDOFF.json','state.rds','seed/READY.json')):
        raise RuntimeError('The earlier online workflow has committed state; inspect it before switching')

def inspect_controller(project, job):
    if not job.isdigit():raise RuntimeError('Expected one scalar controller job ID')
    output=subprocess.check_output(['scontrol','show','job','-o',job],universal_newlines=True,timeout=15)
    fields=dict(re.findall(r'(?:^|\s)(\w+)=(.*?)(?=\s\w+=|$)',output.strip()))
    expected=str(project/'slider_inplace_runtime_2026-10-29/job/resume_slide_prior_em')
    if fields.get('JobId')!=job or fields.get('JobName')!='slide_inplace_control' or \
       fields.get('Command','').split(' ')[0]!=expected or \
       fields.get('UserId','').split('(')[0]!=os.environ['USER']:
        raise RuntimeError('Controller identity does not match this project and user')
    return fields

def check_queue(rows, cfg, mb):
    allowed=('slide_inplace_control','slide_inplace_fit','slide_inplace_validate')
    if any(r['name'] not in allowed for r in mb.legacy_jobs(rows)):
        raise RuntimeError('An unexpected older slider controller/worker is present')
    if any(r['name'].startswith(('slide_mb_','slide_online_')) and r['name'] not in mb.names(cfg) for r in rows):
        raise RuntimeError('Another online workflow has active jobs')

def begin(branch,mb):
    cfg=mb.read(branch/'config.json');project=Path(cfg['source'])
    root=project/'results_slide_prior_em';journal=branch/'handoff_pending.json'
    mb.verify_source(project);boundary(project,mb)
    deadline=time.monotonic()+24*3600
    record=mb.read(journal) if journal.exists() else dict(operations=[])
    while True:
        boundary(project,mb)
        rows=mb.queue();check_queue(rows,cfg,mb)
        controls=[r for r in rows if r['name']=='slide_inplace_control']
        if not controls:break
        # An old controller can briefly coexist with its pending successor.
        # Wait for the running controller to leave, then hold the successor.
        if len(controls)!=1 or controls[0]['state']!='PENDING':
            if time.monotonic()>deadline:raise RuntimeError('No pending legacy controller within 24 hours')
            print('Waiting for a pending controller; existing fitting jobs continue. '+
                  str([(r['job'],r['state'],r['reason']) for r in controls]),flush=True)
            time.sleep(20);continue
        job=controls[0]['job'];fields=inspect_controller(project,job)
        if fields.get('JobState')!='PENDING':
            time.sleep(20);continue
        if (root/'parallel_continuation_job_id.txt').read_text().strip()!=job:
            if time.monotonic()>deadline:raise RuntimeError('Continuation ID did not match pending controller')
            time.sleep(20);continue
        op=dict(job=job,phase='hold_requested',time=time.time(),details=fields)
        record['operations'].append(op);mb.atomic(journal,record)
        subprocess.check_call(['scontrol','hold',job],timeout=15)
        held=inspect_controller(project,job)
        if held.get('JobState')!='PENDING':
            # Holding a job which just started does not cancel/suspend it.
            # Remove our hold and wait for its successor; never scancel it.
            subprocess.check_call(['scontrol','release',job],timeout=15)
            op['phase']='started_before_hold';mb.atomic(journal,record)
            if time.monotonic()>deadline:raise RuntimeError('Controller kept starting before handoff')
            time.sleep(20);continue
        if held.get('Priority')!='0':raise RuntimeError('Pending controller was not held; nothing cancelled')
        op.update(phase='held',held_details=held);mb.atomic(journal,record)
        boundary(project,mb)
        rows=mb.queue();check_queue(rows,cfg,mb)
        if [r['job'] for r in rows if r['name']=='slide_inplace_control']!=[job]:
            raise RuntimeError('Controller set changed while holding the continuation')
        held=inspect_controller(project,job)
        if held.get('JobState')!='PENDING' or held.get('Priority')!='0':
            raise RuntimeError('Controller is no longer held and pending; nothing cancelled')
        op['phase']='cancel_requested';mb.atomic(journal,record)
        subprocess.check_call(['scancel','--state=PENDING','--name=slide_inplace_control',
                               '--user='+os.environ['USER'],job],timeout=15)
        for attempt in range(15):
            if not any(r['job']==job for r in mb.queue()):break
            time.sleep(2)
        else:raise RuntimeError('Held controller still appears in queue; inspect before resuming')
        op['phase']='cancelled';mb.atomic(journal,record)
        print('Replaced held pending controller '+job+'; every fitting/validation job is preserved.',flush=True)
    boundary(project,mb)
    rows=mb.queue();check_queue(rows,cfg,mb)
    if any(r['name']=='slide_inplace_control' for r in rows):raise RuntimeError('Legacy controller reappeared')
    mb.atomic(branch/'HANDOFF.json',dict(source_iteration=14,first_online_iteration=15,
        last_online_iteration=20,last_full_iteration=25,finisher='online controller with unchanged legacy planner',
        operations=record['operations'],time=time.time()))
    print('Handoff secured: the new controller will finish/certify 014, then start online 015.',flush=True)
