"""Recover a positive worker-array acknowledgement; never infer rejection from silence."""
import os, re, time
from pathlib import Path


def accounting(branch, mb, job, task, require_complete=False):
    cfg=mb.read(branch/'config.json')
    command=['sacct','-n','-P','-X','--array','-j',job,
             '--format=JobID%40,JobName%100,User%128,State%40,ExitCode,WorkDir%1024']
    output=mb.queue_output(command,os.environ.copy())
    seen=set()
    for line in output.splitlines():
        if not line.strip():continue
        fields=[x.strip() for x in line.split('|')]
        if len(fields)!=6:raise RuntimeError('Unexpected array accounting response')
        jid,name,user,state,exitcode,directory=fields
        m=re.fullmatch(re.escape(job)+r'_(\d+)',jid)
        if not m or int(m[1]) in seen:
            raise RuntimeError('Unexpected/duplicate accounting job '+jid)
        seen.add(int(m[1]))
        if name!=mb.names(cfg)[1] or user!=os.environ['USER'] or Path(directory).resolve()!=branch.resolve():
            raise RuntimeError('Accepted array identity mismatch: '+line)
        if require_complete and (state!='COMPLETED' or exitcode!='0:0'):
            raise RuntimeError('Accepted array has not completed successfully: '+line)
    if seen!=set(range(1,len(task['genes'])+1)):
        return None  # Accounting may lag; absence is never proof of rejection.
    return output


def recover(branch, mb, pending, expected_job=None, require_complete=False):
    if pending.get('kind')!='workers' or not pending.get('task'):
        raise RuntimeError('Only worker-array acknowledgements can be reconciled here')
    task=pending['task'];td=Path(task['directory'])
    if td.resolve().parent!=(branch/'tasks').resolve() or td.is_symlink() or mb.read(td/'task.json')!=task:
        raise RuntimeError('Pending worker task changed or outside workflow')
    d=Path(task['stage_dir'])
    if mb.sha(d/'input.rds')!=task['input_sha256']:
        raise RuntimeError('Pending worker prior/input changed')
    candidates=set()
    for p in td.glob('worker_*.err'):
        m=re.fullmatch(r'worker_(\d+)_(\d+)\.err',p.name)
        if not m or not 1<=int(m[2])<=len(task['genes']):
            raise RuntimeError('Unrecognized worker log for pending task')
        candidates.add(m[1])
    # A started worker produces a log in this task's unique output directory.
    # Do not guess from an empty queue or resubmit an unacknowledged array.
    if len(candidates)>1:raise RuntimeError('Multiple arrays wrote into the pending task directory')
    if not candidates:return None
    job=candidates.pop()
    if expected_job and job!=expected_job:
        raise RuntimeError('Accepted array differs from the reported job')
    evidence=accounting(branch,mb,job,task,require_complete)
    if evidence is None:return None
    if require_complete:
        for gene in task['genes']:
            if not mb.is_done(branch,d,gene):
                raise RuntimeError('Completed worker lacks a valid saved result: '+gene)
    result=dict(pending,job=job,recovered_acknowledgement=dict(time=time.time(),evidence='unique task worker logs and Slurm accounting'))
    return result,evidence


def wait_for_receipt(branch, mb, pending, seconds=20*3600):
    deadline=time.monotonic()+seconds
    while True:
        recovered=recover(branch,mb,pending)
        if recovered:
            payload,evidence=recovered
            td=Path(payload['task']['directory'])
            (td/'recovered_accounting.txt').write_text(evidence)
            return payload
        if time.monotonic()>=deadline:
            raise RuntimeError('Worker submission still unacknowledged: preserved pending_submission.json; no duplicate array submitted')
        print('Waiting for positive evidence of accepted worker submission; no array is resubmitted.',flush=True)
        time.sleep(60)
