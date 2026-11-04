"""Build isolated maintenance and expanded-slider tools from the deployed bundle."""
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

ROOT = Path(__file__).resolve().parents[2]
SRC = Path(__file__).resolve().parent
NAME = "slider_storage_2026-10-28"
DEST = ROOT / "output" / NAME
DEST.mkdir(parents=True, exist_ok=True)

def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(data.replace("\r\n", "\n"), encoding="utf-8", newline="\n")

with zipfile.ZipFile(ROOT / "output/slide_prior_shared_rcc_2026-10-23.zip") as z:
    originals = {}
    for name in z.namelist():
        if name.startswith("script/scan_tissue_attempt/") and name.endswith(".R") and "/tests/" not in name:
            baseline_text = z.read(name).decode().replace("\r\n", "\n")
            write(DEST / "baseline" / name, baseline_text)
            originals[name] = hashlib.sha256(baseline_text.encode()).hexdigest()
        if name in ("job/resume_slide_prior_em", "job/slide_prior_parallel_worker",
                    "script/scan_tissue_attempt/workhorse_slide_prior_em.R"):
            text = z.read(name).decode().replace("\r\n", "\n")
            originals[name] = hashlib.sha256(text.encode()).hexdigest()
            if name.endswith("workhorse_slide_prior_em.R"):
                old = "    if (!is.null(warm)) args$model_init <- warm"
                assert text.count(old) == 1
                text = text.replace(old, "    if (!is.null(warm)) { warm$call <- NULL; args$model_init <- warm }")
                old = "    fit <- do.call(susieRSlidePrior::susie, args)"
                assert text.count(old) == 1
                text = text.replace(old, old + "\n    # Call metadata embeds genotypes and all previous fits when do.call is used.\n    fit$call <- NULL")
            write(DEST / "patch" / name, text)
    write(DEST / "baseline/script/scan_tissue_attempt/slide_shared_scheduler.py",
          z.read("script/scan_tissue_attempt/slide_shared_scheduler.py").decode())

# Stable-ID workers have already been regression-tested against both models.
old = ROOT / "output/renamed_gene_preflight_2026-10-26"
for f in ("renamed_gene_map.csv", "preflight.R", "check_bim.py"):
    shutil.copyfile(old / f, DEST / f)
for f in ("em_utils.R", "slide_prior_em_utils.R", "workhorse_utils.R", "get_gene_annotations.R", "workhorse_slide_prior_em.R"):
    text = ((DEST / "baseline/script/scan_tissue_attempt" / f) if f in ("em_utils.R", "slide_prior_em_utils.R")
            else (old / "script/scan_tissue_attempt" / f)).read_text()
    if f in ("workhorse_utils.R", "get_gene_annotations.R"):
        originals["script/scan_tissue_attempt/" + f] = hashlib.sha256(text.replace("\r\n", "\n").encode()).hexdigest()
    if f == "workhorse_slide_prior_em.R":
        text = text.replace("    if (!is.null(warm)) args$model_init <- warm",
                            "    if (!is.null(warm)) { warm$call <- NULL; args$model_init <- warm }")
        text = text.replace("    fit <- do.call(susieRSlidePrior::susie, args)",
                            "    fit <- do.call(susieRSlidePrior::susie, args)\n    fit$call <- NULL")
    write(DEST / "script/scan_tissue_attempt" / f, text)
write(DEST / "expected_sources.json", json.dumps(originals, indent=2) + "\n")
for p in SRC.iterdir():
    if p.is_file() and p.name != "build_bundle.py":
        shutil.copyfile(p, DEST / p.name)
for p in (SRC / "tests").glob("*"):
    if p.is_file():
        (DEST / "tests").mkdir(exist_ok=True)
        shutil.copyfile(p, DEST / "tests" / p.name)
files = sorted(p for p in DEST.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
assert not any(set(p.relative_to(DEST).parts) & {"checks", "state", "__pycache__"} for p in files)
write(DEST / "SHA256SUMS", "".join(hashlib.sha256(p.read_bytes()).hexdigest() + "  " + p.relative_to(DEST).as_posix() + "\n" for p in files))
archive = ROOT / "output" / (NAME.replace("storage_", "storage_rcc_") + ".zip")
with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
    for p in files + [DEST / "SHA256SUMS"]:
        z.write(p, NAME + "/" + p.relative_to(DEST).as_posix())
print(archive)
print(hashlib.sha256(archive.read_bytes()).hexdigest())
