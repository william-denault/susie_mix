import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]).resolve()
state=json.loads((root/'frozen_expansion.json').read_text())
for rel,digest in state['hashes'].items():
    if hashlib.sha256((root/rel).read_bytes()).hexdigest()!=digest:
        raise SystemExit('Expanded source or frozen gene map changed: '+rel)
if not (root/'INITIALIZED.json').exists(): raise SystemExit('Expanded initialization is incomplete')
