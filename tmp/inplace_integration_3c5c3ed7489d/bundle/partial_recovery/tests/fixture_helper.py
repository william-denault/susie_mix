import contextlib,hashlib,json,os,sys
from pathlib import Path
sys.dont_write_bytecode=True
bundle=Path(sys.argv[1]).resolve();sys.path.insert(0,str(bundle))
import storage,partial_resume
root=Path(sys.argv[2]).resolve();project=root/'project';original=root/'original_state';scratch=root/'scratch';state=root/'final_state'
if sys.argv[3]=='compact_one':
    scratch.mkdir()
    storage.atomic_json(original/'identity.json',dict(project=str(project),kind='compact',expected=2,scratch=str(scratch)))
    storage.atomic_json(scratch/'owner.json',dict(project=str(project),state=str(original)))
    storage.atomic_json(original/'files.json',storage.inventory(project,2))
    target=project/'script/scan_tissue_attempt/workhorse_slide_prior_em.R';target.parent.mkdir(parents=True)
    target.write_bytes((bundle/'baseline/script/scan_tissue_attempt/workhorse_slide_prior_em.R').read_bytes())
    rel='results_slide_prior_em/iteration_001/results/GENEA.rds'
    stage=scratch/hashlib.sha256(rel.encode()).hexdigest()[:24]
    ledger=storage.Ledger(original/'receipts.sqlite')
    storage.transact(project,rel,stage,ledger);ledger.db.close()
elif sys.argv[3]=='finalize':
    @contextlib.contextmanager
    def local_lock(*a):yield
    partial_resume.exclusive=local_lock
    storage.check_sources=lambda p:None  # Fixture supplies only the worker source.
    partial_resume.finalize(project,original,scratch,state,expected=2)
else:raise RuntimeError('Unknown fixture mode')
