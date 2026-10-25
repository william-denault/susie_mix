"""Build a small atomic-install patch; preserve the running fitter and scheduler."""
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parents[2]
PREFIX = "updates/slide_prior_target_40/"
FILES = ["script/scan_tissue_attempt/slide_prior_run_control.R",
         "script/scan_tissue_attempt/run_slide_prior_parallel.R"]


def normalized(data):
    return data.decode("utf-8-sig").replace("\r\n", "\n").encode("utf-8")


INSTALLER = '''"""Atomically install target controls. Does not submit, cancel or hold jobs."""
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

def digest(data):
    return hashlib.sha256(data.decode("utf-8-sig").replace("\\r\\n", "\\n").encode("utf-8")).hexdigest()

def main():
    if len(sys.argv) != 2:
        raise RuntimeError("Usage: python3 install.py RCC_PROJECT_DIRECTORY")
    project = Path(sys.argv[1]).resolve()
    here = Path(__file__).resolve().parent
    manifest = json.loads((here / "manifest.json").read_text())
    if not (project / "results_slide_prior_em/prior_history.csv").is_file():
        raise RuntimeError("No existing slider EM run at this project path")
    prepared = []
    # Validate every input before replacing any installed file. Order is helper
    # first, then the CLI that loads it; each replacement is a single rename.
    for entry in manifest:
        destination = (project / entry["path"]).resolve()
        if project not in destination.parents:
            raise RuntimeError("Invalid destination")
        content = (here / "payload" / entry["path"]).read_bytes()
        if digest(content) != entry["sha256"]:
            raise RuntimeError("Payload checksum mismatch: " + entry["path"])
        if destination.exists():
            current = destination.read_bytes()
            if digest(current) not in entry["accepted"]:
                raise RuntimeError("Installed file differs from the reviewed version: " + entry["path"])
        elif entry["required"]:
            raise RuntimeError("Required existing CLI is missing")
        else:
            current = None
        prepared.append((destination, content, current))
    for destination, content, current in prepared:
        if current == content:
            continue
        if current is not None:
            backup = here / "backup" / destination.name
            backup.parent.mkdir(parents=True, exist_ok=True)
            if not backup.exists():
                backup.write_bytes(current)
        fd, temporary = tempfile.mkstemp(prefix=".slide_target_", dir=str(destination.parent))
        try:
            with os.fdopen(fd, "wb") as stream:
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
            os.chmod(temporary, 0o644)
            os.replace(temporary, str(destination))
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    print("Installed target controls atomically. No jobs, fits, priors or simulation sources changed.")
    print("Next run the R CLI audit, then set-target PROJECT latest 40.")

if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("STOP: " + str(error), file=sys.stderr)
        sys.exit(1)
'''

with ZipFile(ROOT / "output/slide_prior_shared_rcc_2026-10-23.zip") as old:
    contents = {}
    manifest = []
    for name in FILES:
        data = normalized((ROOT / name).read_bytes())
        sha = hashlib.sha256(data).hexdigest()
        required = name in old.namelist()
        accepted = [sha]
        if required:
            accepted.append(hashlib.sha256(normalized(old.read(name))).hexdigest())
        manifest.append(dict(path=name, sha256=sha, accepted=accepted, required=required))
        contents[PREFIX + "payload/" + name] = data
contents[PREFIX + "install.py"] = INSTALLER.encode("utf-8")
contents[PREFIX + "manifest.json"] = (json.dumps(manifest, indent=2) + "\n").encode("utf-8")
contents[PREFIX + "README.md"] = normalized((ROOT / "job/README_slide_prior_parallel.md").read_bytes())
for name in ("test_slide_prior_run_control.R", "test_slide_prior_em_worker.R"):
    contents[PREFIX + "tests/" + name] = normalized((ROOT / "script/scan_tissue_attempt/tests" / name).read_bytes())
destination = ROOT / "output/slide_prior_target_40_rcc_2026-10-25.zip"
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
