import importlib.util,shutil,sys
from pathlib import Path
sys.dont_write_bytecode=True
src=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(src))
import expansion
project=Path('tmp/slider_storage_fixture_path.txt').read_text().strip()
root=Path(project).parent
bundle=root/'test_bundle'
shutil.copytree('output/slider_storage_2026-10-28',bundle,dirs_exist_ok=True)
branch=root/'expanded'
expansion.BUNDLE=bundle
expansion.make_code(branch)
Path('tmp/slider_expansion_fixture_path.txt').write_text(str(root))
print(root)
