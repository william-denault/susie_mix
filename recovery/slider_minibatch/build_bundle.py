"""Build an independent RCC online-EM experiment; never alter earlier bundles."""
import hashlib, json, shutil, zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / 'output/slider_inplace_646_2026-10-30'
DEST = ROOT / 'output/slider_minibatch_2026-11-04'
MODEL = ['em_utils.R', 'slide_prior_em_utils.R', 'slide_prior_recovery.R',
         'slide_prior_parallel.R', 'workhorse_utils.R', 'get_gene_annotations.R',
         'workhorse_slide_prior_em.R']

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def build():
    DEST.mkdir(parents=True, exist_ok=True)
    for name in MODEL:
        src = BASE / 'script/scan_tissue_attempt' / name
        if not src.exists(): src = BASE / 'baseline/script/scan_tissue_attempt' / name
        dst = DEST / 'script/scan_tissue_attempt' / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(str(src), str(dst))
    for name in ('renamed_gene_map.csv', 'excluded_candidates.csv'):
        shutil.copyfile(str(BASE / name), str(DEST / name))
    for p in HERE.rglob('*'):
        if not p.is_file() or '__pycache__' in p.parts: continue
        d = DEST / p.relative_to(HERE)
        d.parent.mkdir(parents=True, exist_ok=True)
        d.write_bytes(p.read_bytes().replace(b'\r\n', b'\n'))
    files = sorted(p for p in DEST.rglob('*') if p.is_file() and
                   p.name != 'SHA256SUMS' and '__pycache__' not in p.parts)
    # Keep the checksum manifest usable by Linux sha256sum after a Windows build.
    (DEST / 'SHA256SUMS').write_bytes(''.join(digest(p) + '  ' + p.relative_to(DEST).as_posix() + '\n' for p in files).encode('utf-8'))
    archive = ROOT / 'output/slider_minibatch_rcc_2026-11-04.zip'
    with zipfile.ZipFile(str(archive), 'w', zipfile.ZIP_DEFLATED) as z:
        for p in files + [DEST / 'SHA256SUMS']:
            z.write(str(p), DEST.name + '/' + p.relative_to(DEST).as_posix())
    print(str(archive)); print('SHA256: ' + digest(archive))

if __name__ == '__main__': build()
