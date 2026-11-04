"""Finish interrupted file transactions, repair metadata, and resume only 005."""
import argparse,csv,hashlib,json,os,re,subprocess,sys,time
from contextlib import contextmanager
from pathlib import Path
import storage
BUNDLE=Path(__file__).resolve().parent

def unsafe_jobs(rows,current=''):
    return [r for r in rows if r[0]!=current and
        (r[1].startswith(('slide_parallel_','slide_expanded_','slide_recover_','slide_inplace_','slide_storage_','slide_partial_','slide_roll_')) or
         r[1] in ('em_slide_prior_prepare','em_susie_slide_prior','slide_prior_em','slide_expand_setup'))]

@contextmanager
def exclusive(project):
    with storage.maintenance_lock(project):
        live=unsafe_jobs(storage.queue(),os.environ.get('SLURM_JOB_ID',''))
        if live: raise RuntimeError('Maintenance/slider jobs still active; no files changed: '+repr(live[:8]))
        yield

def valid_relative(project,rel):
    if not re.fullmatch(r'results_slide_prior_em/iteration_[0-9]{3,}/results/[A-Za-z0-9_.-]+\.rds',rel):
        raise RuntimeError('Invalid result path: '+rel)
    p=project/rel
    if p.is_symlink() or not storage.contained(p,project) or p.resolve()!=p:
        raise RuntimeError('Redirected result path: '+rel)
    return p

def no_new_compaction(*args):
    raise RuntimeError('Partial finalization must never initiate new compaction')

def recover_transactions(project,scratch,ledger,files):
    recovered=0
    for journal in sorted(scratch.glob('*/transaction.json')):
        t=json.loads(journal.read_text());rel=t['relative']
        if rel not in files or journal.parent.name!=hashlib.sha256(rel.encode()).hexdigest()[:24]:
            raise RuntimeError('Unexpected scratch transaction: '+str(journal))
        src=valid_relative(project,rel)
        if ledger.get(rel) is None and t['stage']=='start':
            # In START the source has never been removed by the original tool.
            # Ignore incomplete scratch output; retain the proven original.
            if not src.is_file() or storage.stamp(src)!=t['original_stamp'] or storage.hashes(src)!=t['old_hash']:
                raise RuntimeError('Unfinished transaction original changed; retain scratch: '+rel)
            t.update(stage='unchanged',new_hash=t['old_hash'],new_bytes=t['old_bytes'],removed_calls=0,
                     partial_finalization='original retained before replacement')
            storage.atomic_json(journal,t)
        # PREPARED/VACATED/REPLACED have already passed scientific round-trip
        # verification. The original transactional code finishes only these.
        storage.transact(project,rel,journal.parent,ledger,compact=no_new_compaction)
        recovered+=1
        print('Resolved in-flight transaction: '+rel,flush=True)
    return recovered

def verify_committed(project,ledger,files):
    records=[]
    for rel,raw in ledger.db.execute('SELECT relative,receipt FROM files'):
        if rel not in files: raise RuntimeError('Receipt outside original inventory')
        receipt=json.loads(raw);p=valid_relative(project,rel)
        if not p.is_file() or storage.stamp(p)!=receipt['stamp']:
            raise RuntimeError('Committed result changed: '+rel)
        for key in ('old_md5','new_md5'):
            if not re.fullmatch('[a-f0-9]{32}',receipt[key]): raise RuntimeError('Invalid committed checksum')
        records.append(receipt)
    return records

def finalize(project,original_state,scratch,state,expected=5):
    project,original_state,scratch,state=[p.resolve() for p in (project,original_state,scratch,state)]
    if state==original_state or storage.contained(state,project) or storage.contained(scratch,project):
        raise RuntimeError('Use a separate HOME finalization state and the original external scratch')
    storage.verify_bundle();storage.check_sources(project)
    identity=dict(project=str(project),kind='compact',expected=expected,scratch=str(scratch))
    if json.loads((original_state/'identity.json').read_text())!=identity:
        raise RuntimeError('Original compaction identity differs')
    if json.loads((scratch/'owner.json').read_text())!=dict(state=str(original_state),project=str(project)):
        raise RuntimeError('Scratch owner differs; use the ORIGINAL scratch path')
    storage.setup_state(project,state,'partial_finalize',dict(original_state=str(original_state),scratch=str(scratch),expected=expected))
    if (state/'COMPLETE.json').exists():
        print('Partial finalization already complete.',flush=True);return
    with exclusive(project):
        for name in ('snapshot.rds','files.json','prior_history.preserved.csv'):
            src,dst=original_state/name,state/name
            if dst.exists():
                if storage.hashes(src)!=storage.hashes(dst): raise RuntimeError('Original maintenance snapshot changed')
            else: storage.copy_verified(src,dst)
        listed=json.loads((state/'files.json').read_text());files=set(listed)
        if len(files)!=len(listed): raise RuntimeError('Duplicate original inventory')
        for rel in files: valid_relative(project,rel)
        ledger=storage.Ledger(original_state/'receipts.sqlite')
        try:
            recovered=recover_transactions(project,scratch,ledger,files)
            records=verify_committed(project,ledger,files)
            if set(storage.inventory(project,expected))!=files:
                raise RuntimeError('Result inventory differs after transaction recovery')
            ledger.export(state/'receipts.csv')
        finally: ledger.db.close()
        print('Keeping %d remaining results unchanged; %d committed transactions; %d in-flight transactions resolved.'%
              (len(files)-len(records),len(records),recovered),flush=True)
        storage.install(project,state)
        subprocess.check_call([os.environ.get('RSCRIPT','Rscript'),'--vanilla',str(BUNDLE/'partial_finalize.R'),
                               str(BUNDLE),str(project),str(state)])
        # This NEW state certifies a finalized mixed compact/original checkpoint.
        # It does not claim that all original results have been compacted.
        storage.atomic_json(state/'COMPLETE.json',dict(project=str(project),latest=expected,files=len(files),
            mode='partial_compaction_finalized',compaction_complete=False,metadata_finalized=True,
            original_state=str(original_state),converted_files=sum(r['removed_calls']>0 for r in records),
            unchanged_files=len(files)-sum(r['removed_calls']>0 for r in records),
            saved_bytes=sum(r['old_bytes']-r['new_bytes'] for r in records),finished=time.time()))
        storage.atomic_json(original_state/'PARTIAL_FINALIZED.json',dict(final_state=str(state),project=str(project),
             instruction='Use the final state to resume 005. Do not restart the original bulk compaction.'))

def resume(project,state):
    project,state=project.resolve(),state.resolve()
    complete=json.loads((state/'COMPLETE.json').read_text())
    if complete.get('project')!=str(project) or not complete.get('metadata_finalized') or complete.get('latest')!=5:
        raise RuntimeError('Partial finalization is not complete for iteration 005')
    if (state/'resume_job.json').exists():
        print('Resume already submitted: '+(state/'resume_job.json').read_text(),flush=True);return
    with exclusive(project):
        marker=state/'resume_submission.pending.json'
        if marker.exists():raise RuntimeError('Ambiguous previous submission; inspect '+str(marker))
        env=os.environ.copy()
        env.update(SUSIE_MIX_PROJECT_DIR=str(project),SUSIE_SLIDE_TARGET_ITERATION='5',SUSIE_SLIDE_WORKERS='298',
          SUSIE_SLIDE_POLL_MINUTES='5',SUSIE_SLIDE_MAX_ITER='5000',SUSIE_SLIDE_RECORD_ELBO='0',
          SUSIE_SLIDE_EXCLUDE_NODES='midway2-0003,midway2-0011')
        dependency=['--dependency=afterany:'+os.environ['SLURM_JOB_ID']] if os.environ.get('SLURM_JOB_ID') else ['--begin=now+1minutes']
        cmd=['sbatch','--parsable']+dependency+['--exclude=midway2-0003,midway2-0011','--chdir='+str(project/'job'),
          '--output='+str(project/'job/slide_parallel_control_%j.out'),'--error='+str(project/'job/slide_parallel_control_%j.err'),
          str(project/'job/resume_slide_prior_em'),'latest']
        storage.atomic_json(marker,dict(command=cmd,target=5))
        job=subprocess.check_output(cmd,env=env,universal_newlines=True).strip().split(';')[0]
        if not job.isdigit():raise RuntimeError('Ambiguous sbatch response; keep submission marker')
        storage.atomic_json(state/'resume_job.json',dict(job=job,target=5));marker.unlink();storage.sync_dir(state)
        print('Submitted controller '+job+' to finish iteration 005 in the existing folder.',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['run','finalize','resume'])
    for name in ('project','original_state','scratch','state'):p.add_argument(name,type=Path)
    a=p.parse_args()
    try:
        if a.mode in ('run','finalize'):finalize(a.project,a.original_state,a.scratch,a.state)
        if a.mode in ('run','resume'):resume(a.project,a.state)
    except Exception as e:
        print('STOP: '+str(e)+'. Preserve BOTH state directories and scratch; no unverified EM was launched.',file=sys.stderr);sys.exit(1)
