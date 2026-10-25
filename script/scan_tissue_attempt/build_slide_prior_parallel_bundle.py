"""Code-only RCC shared scheduling update; never packages results or numbered simulations."""
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

root = Path(__file__).resolve().parents[2]
files = [
    "job/resume_slide_prior_em",
    "job/slide_prior_parallel_worker",
    "job/slide_shared_sim_worker",
    "job/launch_simulation_slide",
    "job/recover_slide_prior_em",
    "job/README_slide_prior_parallel.md",
    "job/README_slide_prior_recovery.md",
    "script/scan_tissue_attempt/em_utils.R",
    "script/scan_tissue_attempt/slide_prior_em_utils.R",
    "script/scan_tissue_attempt/slide_prior_recovery.R",
    "script/scan_tissue_attempt/prepare_slide_prior_em_iteration.R",
    "script/scan_tissue_attempt/workhorse_slide_prior_em.R",
    "script/scan_tissue_attempt/slide_prior_parallel.R",
    "script/scan_tissue_attempt/run_slide_prior_parallel.R",
    "script/scan_tissue_attempt/slide_prior_run_control.R",
    "script/scan_tissue_attempt/slide_shared_scheduler.py",
    "script/sim/write_jobs.R",
    "script/scan_tissue_attempt/tests/test_slide_prior_parallel.R",
    "script/scan_tissue_attempt/tests/test_slide_prior_parallel_smoke.R",
    "script/scan_tissue_attempt/tests/test_slide_prior_parallel_launcher.sh",
    "script/scan_tissue_attempt/tests/test_slide_shared_scheduler.py",
    "script/scan_tissue_attempt/build_slide_prior_parallel_bundle.py",
]
data = {f: (root / f).read_text(encoding="utf-8-sig").replace("\r\n", "\n").encode("utf-8") for f in files}
checksums = "".join(f"{sha256(data[f]).hexdigest()}  {f}\n" for f in files)
manifest = "job/slide_prior_shared_SHA256SUMS"
(root / manifest).write_text(checksums, encoding="utf-8", newline="\n")
data[manifest] = checksums.encode("utf-8")
target = root / "output/slide_prior_shared_rcc_2026-10-23.zip"
target.parent.mkdir(exist_ok=True)
with ZipFile(target, "w", compression=ZIP_DEFLATED) as z:
    for name, content in data.items():
        info = ZipInfo(name)
        info.create_system = 3
        executable = name.endswith(".sh") or (name.startswith("job/") and not Path(name).suffix and name != manifest)
        info.external_attr = (0o100755 if executable else 0o100644) << 16
        info.compress_type = ZIP_DEFLATED
        z.writestr(info, content)
with ZipFile(target) as z:
    assert z.testzip() is None
    assert set(z.namelist()) == set(data)
    for name, content in data.items():
        assert z.read(name) == content and b"\r" not in content
    assert not any(name.startswith(("results", "script/sim/jobs_slide/", "simulation results/")) for name in z.namelist())
print(f"Verified {len(data)} code/document/test files: {target} ({target.stat().st_size:,} bytes)")
print(f"SHA256 {sha256(target.read_bytes()).hexdigest()}")
