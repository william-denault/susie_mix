"""Build the opt-in tissue-worker patch without rebuilding deployed bundles."""
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parents[2]
PREFIX = "updates/slide_prior_tissue_parallel/"
FILES = [
    "script/scan_tissue_attempt/slide_prior_run_control.R",
    "script/scan_tissue_attempt/slide_prior_tissue_parallel.R",
    "script/scan_tissue_attempt/workhorse_slide_prior_em.R",
    "script/scan_tissue_attempt/slide_prior_parallel.R",
    "job/slide_prior_parallel_worker",
    "script/scan_tissue_attempt/run_slide_prior_parallel.R",
]

def normalized(data):
    return data.decode("utf-8-sig").replace("\r\n", "\n").encode("utf-8")

def digest(data):
    return hashlib.sha256(normalized(data)).hexdigest()

INSTALLER = r'''"""Install optional tissue jobs. No job submission, cancellation or activation."""
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

def digest(data):
    return hashlib.sha256(data.decode("utf-8-sig").replace("\r\n", "\n").encode("utf-8")).hexdigest()

def main():
    if len(sys.argv) != 2:
        raise RuntimeError("Usage: python3 install.py RCC_PROJECT_DIRECTORY")
    project = Path(sys.argv[1]).resolve()
    here = Path(__file__).resolve().parent
    manifest = json.loads((here / "manifest.json").read_text())
    if not (project / "results_slide_prior_em/prior_history.csv").is_file():
        raise RuntimeError("No existing slider EM history at this path")
    prepared = []
    for entry in manifest:
        destination = (project / entry["path"]).resolve()
        if project not in destination.parents:
            raise RuntimeError("Invalid destination")
        content = (here / "payload" / entry["path"]).read_bytes()
        if digest(content) != entry["sha256"]:
            raise RuntimeError("Payload checksum mismatch: " + entry["path"])
        current = destination.read_bytes() if destination.exists() else None
        if current is not None and digest(current) not in entry["accepted"]:
            raise RuntimeError("Installed file differs from reviewed versions: " + entry["path"])
        if current is None and entry["required"]:
            raise RuntimeError("Required installed file missing: " + entry["path"])
        prepared.append((entry, destination, content, current))
    enabled = project / "results_slide_prior_em/tissue_split_genes.rds"
    if enabled.exists() and any(current is None or digest(current) != digest(content)
                                for _, _, content, current in prepared):
        raise RuntimeError("Tissue splitting is already enabled; refuse to change its frozen sources")
    for entry, destination, content, current in prepared:
        if current == content:
            continue
        if current is not None:
            backup = here / "backup" / entry["path"]
            backup.parent.mkdir(parents=True, exist_ok=True)
            if not backup.exists():
                backup.write_bytes(current)
        destination.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(prefix=".slide_tissue_", dir=str(destination.parent))
        try:
            with os.fdopen(fd, "wb") as stream:
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
            os.chmod(temporary, 0o755 if entry["path"].startswith("job/") else 0o644)
            os.replace(temporary, str(destination))
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    print("Installed optional tissue workers. No genes enabled; no jobs or results changed.")
    print("The existing iteration target and shared scheduler are unchanged. See README.md for preflight and selection.")

if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("STOP: " + str(error), file=sys.stderr)
        sys.exit(1)
'''

def build():
    contents, manifest = {}, []
    with ZipFile(ROOT / "output/slide_prior_shared_rcc_2026-10-23.zip") as old, \
         ZipFile(ROOT / "output/slide_prior_target_40_rcc_2026-10-25.zip") as target:
        for name in FILES:
            data = normalized((ROOT / name).read_bytes())
            accepted = {digest(data)}
            required = name in old.namelist()
            if required:
                accepted.add(digest(old.read(name)))
            alternative = "updates/slide_prior_target_40/payload/" + name
            if alternative in target.namelist():
                accepted.add(digest(target.read(alternative)))
            manifest.append(dict(path=name, sha256=digest(data), accepted=sorted(accepted), required=required))
            contents[PREFIX + "payload/" + name] = data
    contents[PREFIX + "install.py"] = INSTALLER.encode("utf-8")
    contents[PREFIX + "manifest.json"] = (json.dumps(manifest, indent=2) + "\n").encode("utf-8")
    contents[PREFIX + "README.md"] = normalized((ROOT / "job/README_slide_prior_tissue_parallel.md").read_bytes())
    for name in ("test_slide_prior_tissue_parallel.R", "test_slide_prior_em_worker.R",
                 "test_slide_prior_parallel_launcher.sh"):
        contents[PREFIX + "tests/" + name] = normalized((ROOT / "script/scan_tissue_attempt/tests" / name).read_bytes())
    destination = ROOT / "output/slide_prior_tissue_parallel_rcc_2026-10-26.zip"
    with ZipFile(destination, "w", compression=ZIP_DEFLATED) as archive:
        for name, content in contents.items():
            info = ZipInfo(name)
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            info.compress_type = ZIP_DEFLATED
            archive.writestr(info, content)
    with ZipFile(destination) as archive:
        assert archive.testzip() is None
        assert all(archive.read(name) == data for name, data in contents.items())
    print(destination)
    print("SHA256 " + hashlib.sha256(destination.read_bytes()).hexdigest())

if __name__ == "__main__":
    build()
