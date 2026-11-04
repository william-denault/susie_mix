"""Create a NEW bundle; never modify the already submitted storage bundle."""
import hashlib
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[2]
SRC = Path(__file__).resolve().parent
NAME = 'slider_rolling_2026-10-28'
DEST = ROOT / 'output' / NAME
BASE = ROOT / 'output/slider_storage_rcc_2026-10-28.zip'
EXCLUDED = {'submit_recovery.sh', 'storage_job'}
assert hashlib.sha256(BASE.read_bytes()).hexdigest() == 'cdc503c58ded5244e2167359a47f16a67fe7ba8d8cf37f31634b239a9c6ab9b0'

def write(p, data):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(data.replace('\r\n', '\n').encode())

with zipfile.ZipFile(BASE) as z:
    for member in z.namelist():
        rel = '/'.join(member.split('/')[1:])
        if not rel or rel == 'SHA256SUMS' or rel in EXCLUDED or rel.startswith('tests/'):
            continue
        write(DEST / rel, z.read(member).decode())
for p in SRC.rglob('*'):
    if p.is_file() and p.name != 'build_bundle.py' and '__pycache__' not in p.parts:
        write(DEST / p.relative_to(SRC), p.read_text())
files = sorted(p for p in DEST.rglob('*') if p.is_file() and p.name != 'SHA256SUMS' and p.name not in EXCLUDED)
assert not any(set(p.relative_to(DEST).parts) & {'checks', 'state', '__pycache__'} for p in files)
write(DEST / 'SHA256SUMS', ''.join(hashlib.sha256(p.read_bytes()).hexdigest() + '  ' + p.relative_to(DEST).as_posix() + '\n' for p in files))
archive = ROOT / 'output' / (NAME.replace('rolling_', 'rolling_rcc_') + '.zip')
with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as z:
    for p in files + [DEST / 'SHA256SUMS']:
        z.write(p, NAME + '/' + p.relative_to(DEST).as_posix())
print(str(archive))
print(hashlib.sha256(archive.read_bytes()).hexdigest())
