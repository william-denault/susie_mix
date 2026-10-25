"""Exercise the exact shipped installer against temporary RCC-shaped directories."""
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
from zipfile import ZipFile

repo = Path(__file__).resolve().parents[3]
patch_zip = repo / "output/slide_prior_target_40_rcc_2026-10-25.zip"
cli = "script/scan_tissue_attempt/run_slide_prior_parallel.R"
helper = "script/scan_tissue_attempt/slide_prior_run_control.R"
with tempfile.TemporaryDirectory(prefix="target_patch_", dir=str(repo / "tmp")) as temporary:
    root = Path(temporary).resolve()
    assert (repo / "tmp").resolve() in root.parents
    with ZipFile(repo / "output/slide_prior_shared_rcc_2026-10-23.zip") as old:
        old.extractall(root)
    run = root / "results_slide_prior_em"
    run.mkdir()
    (run / "prior_history.csv").write_text("sentinel: do not modify\n")
    (run / "result.rds").write_bytes(b"existing scientific result")
    before = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in root.rglob("*") if p.is_file() and str(p.relative_to(root)).replace("\\", "/") != cli}
    with ZipFile(patch_zip) as archive:
        archive.extractall(root)
    install = root / "updates/slide_prior_target_40/install.py"
    for _ in range(2):
        subprocess.check_call([sys.executable, str(install), str(root)])
    for name, expected in before.items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
    for name in (cli, helper):
        assert (root / name).read_bytes() == (repo / name).read_text(encoding="utf-8-sig").replace("\r\n", "\n").encode()
    # Any unfamiliar installed edit is rejected before another file is changed.
    (root / cli).write_text("unreviewed customization\n")
    helper_before = (root / helper).read_bytes()
    answer = subprocess.run([sys.executable, str(install), str(root)], universal_newlines=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert answer.returncode == 1 and "differs from the reviewed version" in answer.stderr
    assert (root / helper).read_bytes() == helper_before
    assert (root / cli).read_text() == "unreviewed customization\n"
print("PASS: shipped installer, idempotence, unfamiliar-version rejection, all original data and fitter/scheduler files unchanged.")
