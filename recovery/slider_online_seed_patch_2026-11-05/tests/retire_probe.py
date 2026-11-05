"""Used only by the tiny real-R test, on its temporary workspace fixture."""
import importlib.util,json,sys
from pathlib import Path
from unittest.mock import patch
base=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('probe_mb',base/'runtime/minibatch.py')
mb=importlib.util.module_from_spec(spec);spec.loader.exec_module(mb)
branch=Path(sys.argv[1]);project=branch.parents[1]
if not project.name.startswith('online_real_') or Path('tmp').resolve() not in project.resolve().parents:
    raise RuntimeError('This probe only operates on a local temporary test fixture')
mb.atomic(branch/'config.json',dict(source=str(project),source_iteration=14,id='test'))
plan=mb.csvrows(branch/'stages.csv');row=plan[int(sys.argv[2])-1]
for item in plan:
    d=project/'results_slide_prior_em'/('iteration_%03d'%int(item['iteration']))
    d.mkdir(exist_ok=True);mb.atomic(d/'ONLINE_OWNER.json',dict(workflow_id='test',iteration=int(item['iteration'])))
with patch.object(mb,'queue',return_value=[]):mb.retire(branch,row,plan)
