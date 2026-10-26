"""Check the shipped installer against the deployed bundle and target-only patch."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from zipfile import ZipFile

repo = Path(__file__).resolve().parents[3]
patch = repo / "output/slide_prior_tissue_parallel_rcc_2026-10-26.zip"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run_installer(root, good=True):
    answer = subprocess.run([sys.executable, str(root / "updates/slide_prior_tissue_parallel/install.py"), str(root)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    assert (answer.returncode == 0) == good, answer.stdout + answer.stderr
    return answer

for use_target in (False, True):
    with tempfile.TemporaryDirectory(prefix="tissue_patch_", dir=str(repo / "tmp")) as temporary:
        root = Path(temporary).resolve()
        assert (repo / "tmp").resolve() in root.parents
        with ZipFile(repo / "output/slide_prior_shared_rcc_2026-10-23.zip") as old:
            old.extractall(root)
        if use_target:
            with ZipFile(repo / "output/slide_prior_target_40_rcc_2026-10-25.zip") as target:
                for name in target.namelist():
                    prefix = "updates/slide_prior_target_40/payload/"
                    if name.startswith(prefix):
                        (root / name[len(prefix):]).write_bytes(target.read(name))
        run = root / "results_slide_prior_em"; run.mkdir()
        (run / "prior_history.csv").write_text("sentinel history\n")
        (run / "result.rds").write_bytes(b"scientific result")
        (run / "run_control.rds").write_bytes(b"existing target 20")
        (root / "job/active_worker_sentinel").write_text("do not cancel\n")
        with ZipFile(patch) as archive:
            archive.extractall(root)
        manifest = json.loads((root / "updates/slide_prior_tissue_parallel/manifest.json").read_text())
        changed = {entry["path"] for entry in manifest}
        before = {p.relative_to(root).as_posix(): sha(p) for p in root.rglob("*") if p.is_file()
                  and not p.relative_to(root).as_posix().startswith("updates/")
                  and p.relative_to(root).as_posix() not in changed}
        run_installer(root); run_installer(root)
        for name, expected in before.items():
            assert sha(root / name) == expected, name
        for entry in manifest:
            name = entry["path"]
            assert (root / name).read_bytes() == (repo / name).read_text(encoding="utf-8-sig").replace("\r\n", "\n").encode()
        assert not (run / "tissue_split_genes.rds").exists()
        # Unknown edits must fail before another file is written.
        cli = root / "script/scan_tissue_attempt/run_slide_prior_parallel.R"
        cli.write_text("unreviewed custom source\n")
        helper = root / "script/scan_tissue_attempt/slide_prior_tissue_parallel.R"
        helper_before = sha(helper)
        answer = run_installer(root, good=False)
        assert "differs from reviewed" in answer.stderr and sha(helper) == helper_before
print("PASS: exact shipped installer on both deployed versions, idempotence, no activation or target changes, unfamiliar-source rejection and unrelated files preserved.")
