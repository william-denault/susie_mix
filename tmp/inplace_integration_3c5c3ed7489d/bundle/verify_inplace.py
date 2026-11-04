import hashlib,json,sys
from pathlib import Path
from common import RUNTIME,META
project=Path(sys.argv[1]).resolve()
frozen=json.loads((project/RUNTIME/'frozen.json').read_text())
if frozen['project']!=str(project):raise SystemExit('Wrong in-place runtime project')
for rel,digest in frozen['hashes'].items():
    if hashlib.sha256((project/rel).read_bytes()).hexdigest()!=digest:
        raise SystemExit('Frozen in-place code or input changed: '+rel)
if not (project/META/'INITIALIZED.json').exists():raise SystemExit('In-place setup is unfinished')
