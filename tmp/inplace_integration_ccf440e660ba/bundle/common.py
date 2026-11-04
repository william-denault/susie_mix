import json,os
from contextlib import contextmanager
from pathlib import Path
import storage
RUNTIME='slider_inplace_runtime_2026-10-29'
META='results_slide_prior_em/inplace_expansion_2026-10-29'
CLEANUP='results_slide_prior_em/rolling_cleanup_inplace'

def unsafe(rows,current=''):
    return [r for r in rows if r[0]!=current and
      (r[1].startswith(('slide_parallel_','slide_expanded_','slide_inplace_','slide_recover_','slide_storage_','slide_partial_','slide_roll_')) or
       r[1] in ('em_slide_prior_prepare','em_susie_slide_prior','slide_prior_em','slide_expand_setup','slide_next_setup'))]

@contextmanager
def exclusive(project):
    with storage.maintenance_lock(project):
        live=unsafe(storage.queue(),os.environ.get('SLURM_JOB_ID',''))
        if live:raise RuntimeError('Slider/maintenance jobs still active; no changes: '+repr(live[:8]))
        yield

def recovery_ready(project,state):
    x=json.loads((state/'COMPLETE.json').read_text())
    if x.get('project')!=str(project) or x.get('latest')!=5 or 'saved_bytes' not in x:
        raise RuntimeError('Recovery state is not finalized for this project/005')
    if x.get('mode')=='partial_compaction_finalized' and not x.get('metadata_finalized'):
        raise RuntimeError('Partial-compaction metadata is unfinished')

def check_dependencies(project):
    p=project/'results_slide_prior_em/storage_dependencies.json'
    if p.exists() and json.loads(p.read_text()):
        raise RuntimeError('Another branch still depends on these results; rolling cleanup cannot proceed')
