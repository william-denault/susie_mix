"""Export the current self-contained code, without relying on an older bundle."""
import hashlib, sys, zipfile
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
DEST=ROOT/'output/slider_online_current_code'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
    if DEST.exists() or HERE==DEST.resolve() or HERE in DEST.resolve().parents:
        raise RuntimeError('Choose a new destination outside the source code directory')
    archive=DEST.with_name(DEST.name+'.zip')
    if archive.exists():raise RuntimeError('Destination ZIP already exists')
    DEST.mkdir(parents=True)
    for p in HERE.rglob('*'):
        if p.is_file() and '__pycache__' not in p.parts and p!=HERE/'SHA256SUMS':
            d=DEST/p.relative_to(HERE);d.parent.mkdir(parents=True,exist_ok=True)
            d.write_bytes(p.read_bytes().replace(b'\r\n',b'\n'))
    files=sorted(p for p in DEST.rglob('*') if p.is_file() and p.name!='SHA256SUMS' and '__pycache__' not in p.parts)
    # Linux sha256sum treats a CR before LF as part of the filename. Write
    # bytes so building on Windows cannot translate these newlines to CRLF.
    (DEST/'SHA256SUMS').write_bytes(''.join(sha(p)+'  '+p.relative_to(DEST).as_posix()+'\n' for p in files).encode('utf-8'))
    with zipfile.ZipFile(str(archive),'w',zipfile.ZIP_DEFLATED) as z:
        for p in files+[DEST/'SHA256SUMS']:z.write(str(p),DEST.name+'/'+p.relative_to(DEST).as_posix())
    print(str(archive));print('SHA256: '+sha(archive))
if __name__=='__main__':
    if len(sys.argv)!=2:raise SystemExit('Usage: python3 build_bundle.py NEW_DESTINATION_FOLDER')
    DEST=Path(sys.argv[1]).resolve()
    build()
