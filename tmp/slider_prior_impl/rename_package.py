from pathlib import Path

root = Path(r"C:\Document\Serieux\Travail\Package\git\susieR")
old, new = b"susieSlide", b"susieRSlidePrior"
assert b"Package: " + old in (root / "DESCRIPTION").read_bytes()
files = [root / name for name in ("DESCRIPTION", "NAMESPACE", "README.md", "Makefile")]
for directory in ("R", "man", "tests", "inst", "vignettes"):
    files.extend(p for p in (root / directory).rglob("*")
                 if p.is_file() and p.suffix in (".R", ".Rd", ".Rmd", ".md"))
files.append(root / "validation/audit-root.R")
benchmark = root / "inst/examples/benchmark_slider_prior.R"
changed = []
for path in files:
    if path == benchmark:
        continue
    before = path.read_bytes()
    after = before.replace(old, new)
    if after != before:
        path.write_bytes(after)
        changed.append(str(path.relative_to(root)))

old_help = root / "man/susieSlide-package.Rd"
new_help = root / "man/susieRSlidePrior-package.Rd"
assert not new_help.exists()
old_help.rename(new_help)

ignore = root / ".gitignore"
before = ignore.read_bytes()
if b"susieRSlidePrior.Rcheck" not in before:
    ignore.write_bytes(before.replace(b"susieSlide.Rcheck", b"susieSlide.Rcheck\nsusieRSlidePrior.Rcheck"))

before = benchmark.read_bytes()
after = before.replace(
    b"# susieSlide. Workers are sequential:",
    b"# susieRSlidePrior. Workers are sequential:")
after = after.replace(b"library containing susieSlide >= 0.3.0", b"library containing susieRSlidePrior >= 0.3.0")
after = after.replace(
    b'ns <- loadNamespace("susieSlide",lib.loc=if(method=="susie_slide") legacy_lib else prior_lib)',
    b'ns <- loadNamespace(if(method=="susie_slide") "susieSlide" else "susieRSlidePrior",\n'
    b'                        lib.loc=if(method=="susie_slide") legacy_lib else prior_lib)')
assert after != before
benchmark.write_bytes(after)
changed.append(str(benchmark.relative_to(root)))

readme = root / "README.md"
before = readme.read_bytes()
after = before.replace(
    b"sliders. Ordinary susieR can be installed alongside this package for comparison.",
    b"sliders. Both `susieR` and the original `susieSlide` can be installed alongside\n"
    b"`susieRSlidePrior` for comparison.")
assert after != before
readme.write_bytes(after)
print("Updated package identity in", len(changed), "files and renamed package help.")
print("\n".join(changed))
