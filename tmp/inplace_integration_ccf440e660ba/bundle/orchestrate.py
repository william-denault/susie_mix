"""Start/reuse the partial recovery, then queue the same-folder continuation."""
import argparse,json,subprocess,sys
from pathlib import Path
import storage
BUNDLE=Path(__file__).resolve().parent

def recovery_choice(original,partial):
    if (partial/'handoff.pending.json').exists():raise RuntimeError('Resolve ambiguous partial-recovery handoff before continuing')
    if (partial/'handoff.json').exists():
        if json.loads((partial/'handoff.json').read_text()).get('phase')!='released':
            raise RuntimeError('Partial finalizer is still held or handoff is incomplete; inspect handoff.json')
        return partial,False
    if (partial/'COMPLETE.json').exists():return partial,False
    if (original/'COMPLETE.json').exists():return original,False
    return partial,True

def run(a):
    storage.verify_bundle()
    original=a.original_state.resolve();partial=a.partial_state.resolve()
    recovery,need_handoff=recovery_choice(original,partial)
    if need_handoff:
        subprocess.check_call([sys.executable,str(BUNDLE/'partial_recovery/handoff.py'),str(a.project.resolve()),
          str(original),str(a.scratch.resolve()),str(partial)])
    cmd=[sys.executable,str(BUNDLE/'launch.py'),'submit',str(a.project.resolve()),str(recovery),str(a.launch_state.resolve())]
    if a.resume:cmd.append('--resume')
    subprocess.check_call(cmd)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('project','original_state','scratch','partial_state','launch_state'):p.add_argument(name,type=Path)
    p.add_argument('--resume',action='store_true');a=p.parse_args()
    try:run(a)
    except Exception as e:print('STOP: '+str(e)+'. Inspect existing receipts before retrying.',file=sys.stderr);sys.exit(1)
