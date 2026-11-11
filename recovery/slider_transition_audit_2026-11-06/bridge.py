"""Exercise actual Python cleanup against R test fixtures, never RCC data."""
import argparse, hashlib, importlib.util, json, os, sys
from pathlib import Path
from unittest.mock import patch

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
SCRATCH=(ROOT/'tmp').resolve()

def load(code):
    spec=importlib.util.spec_from_file_location('transition_runtime',code/'minibatch.py')
    mb=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mb;spec.loader.exec_module(mb)
    return mb

def paths(branch,mb):
    rows=mb.csvrows(branch/'audit_warm.csv')
    return {Path(r[k]).resolve() for r in rows for k in ('path','priors')}

def metadata(root,branch):
    # Every retained non-result analysis file, including histories and counts.
    files={p for p in root.rglob('*') if p.is_file() and 'results' not in p.relative_to(root).parts}
    files.update(branch/n for n in ('state.rds','progress.csv','stages.csv','model.rds','original_prior_history.csv'))
    return {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}

def run(code,branch,mode,index):
    branch=branch.resolve();code=code.resolve()
    if SCRATCH not in branch.parents or not branch.parent.parent.name.startswith('transition_'):
        raise RuntimeError('This audit only operates on its own tmp/transition_* fixture')
    project=branch.parent.parent;root=project/'results_slide_prior_em'
    mb=load(code)
    if mode=='init':
        n=len(mb.csvrows(branch/'seed/snapshot_001/expected.csv'))
        tissues=len(mb.csvrows(branch/'seed/snapshot_001/source_priors.csv'))
        mb.atomic(branch/'config.json',dict(source=str(project),source_iteration=14,id='fixture',workers=298,seed=20261104,audit_gene_count=n))
        mb.atomic(branch/'seed/READY.json',dict(snapshot=str(branch/'seed/snapshot_001'),source_iteration=14,genes=n,tissue_fits=n*tissues))
        mb.atomic(branch/'frozen_inputs.json',{})
        return
    stages=mb.csvrows(branch/'stages.csv');row=stages[index-1]
    if mode=='own':
        mb.assert_owned_iteration(branch,int(row['iteration']));return
    if mode=='seed':
        with patch.object(mb,'queue',return_value=[]):
            assert mb.seed_run(branch)
        return
    before=metadata(root,branch);warm=paths(branch,mb)
    warm_hashes={str(p):mb.sha(p) for p in warm}
    if mode=='partial':
        actual_unlink=Path.unlink;count=[0]
        def interrupt(p,*a,**kw):
            actual_unlink(p,*a,**kw)
            if p.suffix=='.rds':
                count[0]+=1
                if count[0]==3:raise KeyboardInterrupt('simulated controller termination during cleanup')
        try:
            with patch.object(mb,'queue',return_value=[]),patch.object(Path,'unlink',interrupt):
                mb.retire(branch,row,stages)
        except KeyboardInterrupt:
            assert count[0]==3
        else:
            raise AssertionError('Cleanup interruption did not occur')
        assert {str(p):mb.sha(p) for p in warm}==warm_hashes
        # The next controller resumes the same cleanup safely.
    with patch.object(mb,'queue',return_value=[]):
        mb.retire(branch,row,stages)
        mb.retire(branch,row,stages)
    assert {str(p):mb.sha(p) for p in warm}==warm_hashes
    after=metadata(root,branch)
    assert all(after.get(p)==digest for p,digest in before.items()), 'Cleanup changed history/counts/state'
    iteration=int(row['iteration'])
    for i in range(iteration-1):
        d=root/('iteration_%03d'%i)/'results'
        assert not list(d.glob('*.rds')), 'Old full results remain: '+str(d)
    n=mb.read(branch/'config.json')['audit_gene_count']
    assert len(list((root/('iteration_%03d'%(iteration-1))/'results').glob('*.rds')))==n
    if iteration==15:
        assert len(list((branch/'seed/snapshot_001/results').glob('*.rds')))==n
    else:
        assert not list((branch/'seed/snapshot_001/results').glob('*.rds'))
    print('PASS cleanup at stage %d / iteration %03d: warm starts + histories intact; current/previous results retained.'%(index,iteration))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('code',type=Path);p.add_argument('branch',type=Path)
    p.add_argument('mode',choices=['init','own','cleanup','partial','seed']);p.add_argument('stage',type=int)
    a=p.parse_args();run(a.code,a.branch,a.mode,a.stage)
