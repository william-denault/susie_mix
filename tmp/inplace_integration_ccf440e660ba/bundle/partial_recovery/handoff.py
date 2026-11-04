"""Replace only the recorded compaction and pending resume with a finalizer."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
import storage
from partial_resume import unsafe_jobs
BUNDLE=Path(__file__).resolve().parent

def identify(original_state,rows):
    entries={}
    for line in (original_state/'queued_jobs.tsv').read_text().splitlines():
        fields=line.split('\t')
        if len(fields)!=2 or not fields[1].isdigit() or fields[0] in entries: raise RuntimeError('Invalid original submission receipt')
        entries[fields[0]]=fields[1]
    if not all(k in entries for k in ('compact','resume')):raise RuntimeError('Original compact/resume job IDs missing')
    active={r[0]:r for r in rows}
    stop=[]
    for key,name in (('resume','slide_storage_resume'),('compact','slide_storage_compact')):
        job=entries[key]
        if job in active:
            row=active[job]
            if row[1]!=name:raise RuntimeError('Job identity differs: '+repr(row))
            if key=='resume' and row[2]!='PENDING':raise RuntimeError('The old resume is already active; inspect its controller before switching')
            stop.append(job)
    competing=unsafe_jobs([r for r in rows if r[0] not in stop])
    if competing:raise RuntimeError('Other slider/maintenance jobs active; no cancellation: '+repr(competing[:8]))
    return entries,stop

def run(a):
    import fcntl
    project,original,scratch,state=[p.resolve() for p in (a.project,a.original_state,a.scratch,a.state)]
    storage.verify_bundle();storage.check_sources(project)
    # Read-only checks happen before the held job or any cancellation.
    identity=json.loads((original/'identity.json').read_text())
    if identity!=dict(project=str(project),kind='compact',expected=5,scratch=str(scratch)):
        raise RuntimeError('Wrong original state/project/scratch')
    if json.loads((scratch/'owner.json').read_text())!=dict(project=str(project),state=str(original)):
        raise RuntimeError('Wrong scratch owner')
    if any(not (original/name).is_file() for name in ('snapshot.rds','files.json','receipts.sqlite')):
        raise RuntimeError('Original compaction snapshot or receipt database is missing')
    if (original/'COMPLETE.json').exists():raise RuntimeError('Original compaction already finished; inspect the queued resume instead')
    if storage.contained(state,project) or state==original:raise RuntimeError('Use a separate HOME finalization state')
    state.mkdir(parents=True,exist_ok=True);(state/'logs').mkdir(exist_ok=True)
    with (state/'.handoff.lock').open('a') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        if (state/'handoff.json').exists() or (state/'handoff.pending.json').exists():
            raise RuntimeError('Handoff receipt exists; inspect it instead of submitting twice')
        entries,stop=identify(original,storage.queue())
        cmd=['sbatch','--parsable','--hold','--exclude=midway2-0003,midway2-0011',
             '--output='+str(state/'logs/finalize_%j.out'),'--error='+str(state/'logs/finalize_%j.err')]
        if stop:cmd+=['--dependency=afterany:'+':'.join(stop)]
        cmd += [str(BUNDLE/'partial_job'),str(BUNDLE),str(project),str(original),str(scratch),str(state)]
        storage.atomic_json(state/'handoff.pending.json',dict(command=cmd,old_jobs=entries))
        job=subprocess.check_output(cmd,universal_newlines=True).strip().split(';')[0]
        if not job.isdigit():raise RuntimeError('Ambiguous finalizer submission; inspect pending marker')
        receipt=dict(job=job,old_jobs=entries,cancel_requested=stop,phase='held',time=time.time())
        storage.atomic_json(state/'handoff.json',receipt)
        (state/'handoff.pending.json').unlink()
        # Re-read the queue: refuse to cancel a resume that started meanwhile.
        _,still_active=identify(original,[r for r in storage.queue() if r[0]!=job])
        if any(j not in stop for j in still_active):raise RuntimeError('Queue changed; finalizer remains held')
        for old in still_active:
            cmd=['scancel']+(['--state=PENDING'] if old==entries['resume'] else [])+[old]
            subprocess.check_call(cmd)
        receipt['phase']='old_jobs_cancelled';storage.atomic_json(state/'handoff.json',receipt)
        subprocess.check_call(['scontrol','release',job])
        receipt['phase']='released';storage.atomic_json(state/'handoff.json',receipt)
        print('Finalizer '+job+' queued. Only recorded compaction/resume jobs were stopped: '+','.join(still_active))
        print('It will repair partial-compaction metadata, then resume 005 automatically. Logs: '+str(state/'logs'))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('project','original_state','scratch','state'):p.add_argument(name,type=Path)
    a=p.parse_args()
    try:run(a)
    except Exception as e:
        print('STOP: '+str(e)+'. If a finalizer was submitted, it remains held until the handoff is resolved.',file=sys.stderr);sys.exit(1)
