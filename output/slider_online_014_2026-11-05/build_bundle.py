"""Build the 013 -> online 014..019 -> full 020..024 handoff bundle."""
import hashlib, shutil, zipfile
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
BASE=ROOT/'output/slider_minibatch_2026-11-04'
DEST=ROOT/'output/slider_online_014_2026-11-05'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
    DEST.mkdir(exist_ok=True)
    for p in (BASE/'script').rglob('*.R'):
        d=DEST/p.relative_to(BASE);d.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(p),str(d))
    for name in ('renamed_gene_map.csv','excluded_candidates.csv'):shutil.copyfile(str(BASE/name),str(DEST/name))
    for p in HERE.rglob('*'):
        if p.is_file() and '__pycache__' not in p.parts:
            d=DEST/p.relative_to(HERE);d.parent.mkdir(parents=True,exist_ok=True)
            d.write_bytes(p.read_bytes().replace(b'\r\n',b'\n'))
    files=sorted(p for p in DEST.rglob('*') if p.is_file() and p.name!='SHA256SUMS' and '__pycache__' not in p.parts)
    (DEST/'SHA256SUMS').write_text(''.join(sha(p)+'  '+p.relative_to(DEST).as_posix()+'\n' for p in files))
    archive=ROOT/'output/slider_online_014_rcc_2026-11-05.zip'
    with zipfile.ZipFile(str(archive),'w',zipfile.ZIP_DEFLATED) as z:
        for p in files+[DEST/'SHA256SUMS']:z.write(str(p),DEST.name+'/'+p.relative_to(DEST).as_posix())
    print(str(archive));print('SHA256: '+sha(archive))
if __name__=='__main__':build()
