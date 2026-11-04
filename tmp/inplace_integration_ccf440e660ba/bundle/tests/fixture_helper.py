"""Local synthetic tests only; Slurm guards are mocked explicitly."""
import contextlib,sys
from pathlib import Path
sys.dont_write_bytecode=True
bundle=Path(sys.argv[1]).resolve();sys.path.insert(0,str(bundle))
import inplace,retention,storage
from common import META,CLEANUP
mode=sys.argv[2];project=Path(sys.argv[3]).resolve()
if mode=='code':
    (project/META).mkdir(parents=True)
    storage.atomic_json(project/META/'identity.json',dict(project=str(project),fixture=True))
    inplace.make_code(project)
elif mode=='freeze':
    inplace.freeze(project);storage.atomic_json(project/META/'INITIALIZED.json',dict(fixture=True))
elif mode in ('interrupt','step'):
    current=int(sys.argv[4]);state=project/CLEANUP/('before_%03d'%current)
    if mode=='step':
        @contextlib.contextmanager
        def local_lock(*args):yield
        retention.barrier_lock=local_lock;retention.check_queue=lambda:None
        retention.step(project,current,project/'queue.txt')
    else:
        if not (state/'plan.json').exists():retention.plan_files(project,state,current)
        def interrupt(p):raise RuntimeError('simulated interruption after unlink')
        try:retention.execute_plan(project,state,current,after_unlink=interrupt)
        except RuntimeError as e:assert str(e)=='simulated interruption after unlink'
        else:raise AssertionError('No interruption')
else:raise RuntimeError('Unknown mode')
