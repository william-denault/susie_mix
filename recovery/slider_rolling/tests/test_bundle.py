"""Local archive and generated-code checks; run from repository root."""
import ast,hashlib,subprocess,sys,tempfile,zipfile
from pathlib import Path
sys.dont_write_bytecode=True
root=Path('output/slider_rolling_2026-10-28').resolve()
archive=root.parent/'slider_rolling_rcc_2026-10-28.zip'
with zipfile.ZipFile(archive) as z:
    manifest=z.read(root.name+'/SHA256SUMS').decode().splitlines()
    for line in manifest:
        digest,rel=line.split('  ',1)
        data=z.read(root.name+'/'+rel)
        assert hashlib.sha256(data).hexdigest()==digest,rel
        assert (root/rel).read_bytes()==data,rel
        if rel.endswith('.py'): ast.parse(data,filename=rel,feature_version=(3,6))
    assert not any('/submit_recovery.sh' in n or n.endswith('/storage_job') for n in z.namelist())
sys.path.insert(0,str(root));import expansion
branch=Path(tempfile.mkdtemp(prefix='rolling_code_',dir='tmp')).resolve();expansion.make_code(branch)
shell='C:/Users/willi/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/git/usr/bin/sh.exe'
scripts=[root/'expansion_job',root/'gate_job',root/'submit_expanded_040.sh']+list((branch/'job').iterdir())
for p in scripts:
    assert b'\r' not in p.read_bytes(),p
    subprocess.run([shell,'-n',str(p)],check=True)
for name in ('em_utils.R','slide_prior_em_utils.R','prepare_slide_prior_em_iteration.R'):
    assert (branch/'script/scan_tissue_attempt'/name).read_bytes()==(root/'baseline/script/scan_tissue_attempt'/name).read_bytes()
print('PASS: archive hashes, Python 3.6 syntax, shell syntax/LF, generated code, unchanged numerical EM and preparation helpers.')
print('Payload files:',len(manifest),'ZIP bytes:',archive.stat().st_size)
print('SHA256:',hashlib.sha256(archive.read_bytes()).hexdigest())
