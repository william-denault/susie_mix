"""Build a code-only RCC bundle with Linux line endings; no fitted results."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED, ZipInfo

root = Path(__file__).resolve().parents[1]
files = [
    "job/em_susie_slide_prior",
    "job/em_susie_slide_prior_array",
    "job/README_slide_prior_em.md",
    "script/scan_tissue_attempt/em_utils.R",
    "script/scan_tissue_attempt/workhorse_utils.R",
    "script/scan_tissue_attempt/get_gene_annotations.R",
    "script/scan_tissue_attempt/slide_prior_em_utils.R",
    "script/scan_tissue_attempt/prepare_slide_prior_em_iteration.R",
    "script/scan_tissue_attempt/run_slide_prior_em_chunk.R",
    "script/scan_tissue_attempt/workhorse_slide_prior_em.R",
    "script/scan_tissue_attempt/tests/test_slide_prior_em.R",
    "script/scan_tissue_attempt/tests/test_slide_prior_em_smoke.R",
    "script/scan_tissue_attempt/tests/test_slide_prior_em_worker.R",
    "script/scan_tissue_attempt/tests/test_slide_prior_em_launcher.sh",
]
target = root / "output/slide_prior_em_rcc.zip"
target.parent.mkdir(parents=True, exist_ok=True)
with ZipFile(target, "w", compression=ZIP_DEFLATED) as bundle:
    for name in files:
        data = (root / name).read_text(encoding="utf-8-sig").replace("\r\n", "\n").encode("utf-8")
        info = ZipInfo(name)
        info.create_system = 3
        executable = name.endswith(".sh") or (name.startswith("job/") and not Path(name).suffix)
        info.external_attr = (0o100755 if executable else 0o100644) << 16
        info.compress_type = ZIP_DEFLATED
        bundle.writestr(info, data)
with ZipFile(target) as bundle:
    assert bundle.namelist() == files
    assert bundle.testzip() is None
    for name in files:
        assert b"\r" not in bundle.read(name)
        assert bundle.read(name).decode("utf-8") == (root / name).read_text(encoding="utf-8-sig")
print(f"Verified {len(files)} code/document/test files: {target} ({target.stat().st_size:,} bytes)")
