import hashlib
from pathlib import Path
import zipfile
ROOT=Path(__file__).resolve().parents[2];SRC=Path(__file__).resolve().parent
NAME='slider_inplace_2026-10-29';DEST=ROOT/'output'/NAME
BASE=ROOT/'output/slider_storage_rcc_2026-10-28.zip'
assert hashlib.sha256(BASE.read_bytes()).hexdigest()=='cdc503c58ded5244e2167359a47f16a67fe7ba8d8cf37f31634b239a9c6ab9b0'
def write(p,text):p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(text.replace('\r\n','\n').encode())
with zipfile.ZipFile(BASE) as z:
    for member in z.namelist():
        rel='/'.join(member.split('/')[1:])
        if rel.startswith(('baseline/','patch/','script/')) or rel in ('storage.py','maintenance.R','expected_sources.json','renamed_gene_map.csv','preflight.R','check_bim.py'):
            write(DEST/rel,z.read(member).decode())
partial=ROOT/'output/slider_partial_rcc_2026-10-29.zip'
assert hashlib.sha256(partial.read_bytes()).hexdigest()=='8a66817d458153b1bbe40d9585eb2204cab0fca9918b7b8006ab418cd46953ea'
with zipfile.ZipFile(partial) as z:
    for member in z.namelist():
        rel='/'.join(member.split('/')[1:])
        if rel:write(DEST/'partial_recovery'/rel,z.read(member).decode())
for p in SRC.rglob('*'):
    if p.is_file() and p.name!='build_bundle.py' and '__pycache__' not in p.parts:write(DEST/p.relative_to(SRC),p.read_text())
files=sorted(p for p in DEST.rglob('*') if p.is_file() and p!=DEST/'SHA256SUMS')
assert not any(set(p.relative_to(DEST).parts)&{'checks','state','__pycache__'} for p in files)
write(DEST/'SHA256SUMS',''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(DEST).as_posix()+'\n' for p in files))
archive=ROOT/'output/slider_inplace_rcc_2026-10-29.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for p in files+[DEST/'SHA256SUMS']:z.write(p,NAME+'/'+p.relative_to(DEST).as_posix())
print(archive);print(hashlib.sha256(archive.read_bytes()).hexdigest())
