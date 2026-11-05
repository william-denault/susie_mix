import hashlib,zipfile
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
DEST=ROOT/'output/slider_online_submit_patch_2026-11-05'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def manifest(folder):
    files=sorted(p for p in folder.rglob('*') if p.is_file() and p!=folder/'SHA256SUMS' and '__pycache__' not in p.parts)
    (folder/'SHA256SUMS').write_bytes(''.join(sha(p)+'  '+p.relative_to(folder).as_posix()+'\n' for p in files).encode())
def build():
    manifest(HERE/'runtime')
    for p in HERE.rglob('*'):
        if p.is_file() and '__pycache__' not in p.parts and p!=HERE/'SHA256SUMS':
            q=DEST/p.relative_to(HERE);q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes(p.read_bytes().replace(b'\r\n',b'\n'))
    manifest(DEST/'runtime');manifest(DEST)
    archive=ROOT/'output/slider_online_submit_patch_rcc_2026-11-05.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(DEST.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts:z.write(p,DEST.name+'/'+p.relative_to(DEST).as_posix())
    print(archive);print('SHA256: '+sha(archive))
if __name__=='__main__':build()
