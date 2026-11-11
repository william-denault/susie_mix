"""Reconcile the specifically reported, completed array 50014103 in iteration 023."""
import subprocess, sys
from pathlib import Path
from support import read, encoded, digest
sys.path.insert(0,str(Path(__file__).resolve().parent/'runtime'))
import submission_receipts

ARRAY='50014103'
TASK='task_001411'


def prepare(branch,mb):
    errorfile=branch/'ERROR.json';pendingfile=branch/'pending_submission.json'
    if not errorfile.exists() and not pendingfile.exists():return None
    if not errorfile.exists() or not pendingfile.exists():
        raise RuntimeError('Incomplete failure records; refusing to infer the missing submission')
    error=read(errorfile);pending=read(pendingfile);task=pending.get('task',{})
    td=branch/'tasks'/TASK;cfg=read(branch/'config.json')
    d=Path(cfg['source'])/'results_slide_prior_em/iteration_023/batches/full_step'
    if pending.get('kind')!='workers' or 'job' in pending or task.get('stage')!=93 or \
       Path(task.get('directory','')).resolve()!=td.resolve() or Path(task.get('stage_dir','')).resolve()!=d.resolve() or \
       len(task.get('genes',[]))!=130 or len(set(task['genes']))!=130 or task!=read(td/'task.json'):
        raise RuntimeError('Not the reported iteration-023 submission failure')
    command=['sbatch','--parsable','--export=ALL','--time=23:00:00','--exclude=midway2-0003,midway2-0011',
        '--chdir='+str(branch),'--job-name='+mb.names(cfg)[1],
        '--output='+str(td/'worker_%A_%a.out'),'--error='+str(td/'worker_%A_%a.err'),'--array=1-130',
        str(branch/'code/worker_job'),str(branch),str(td/'task.json')]
    if pending['command']!=command or error.get('job')!='50013914' or error.get('error')!=str(subprocess.CalledProcessError(1,command)):
        raise RuntimeError('Different controller failure; no error marker cleared')
    if 'Batch job submission failed: Unexpected message received' not in (branch/'logs/control_50013914.err').read_text(errors='replace'):
        raise RuntimeError('Original submission error log does not match')
    if (td/'submission.json').exists() or (d/'commit.rds').exists():
        raise RuntimeError('Reported task/update already changed')
    recovered=submission_receipts.recover(branch,mb,pending,expected_job=ARRAY,require_complete=True)
    if not recovered:raise RuntimeError('Accounting/log evidence not yet complete for array '+ARRAY)
    payload,evidence=recovered
    return dict(array=ARRAY,workers=130,submission_path='tasks/'+TASK+'/submission.json',
        submission=payload,evidence=evidence,
        archives={name:dict(sha256=digest(branch/name),text=(branch/name).read_text()) for name in ('ERROR.json','pending_submission.json')})
