import importlib.util, json, os, shutil, tempfile, unittest
from types import SimpleNamespace
from pathlib import Path
from unittest.mock import patch
from contextlib import contextmanager

HERE=Path(__file__).resolve().parents[1]
def load(name):
    spec=importlib.util.spec_from_file_location(name,HERE/(name+'.py'))
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
bim=load('check_bim');backup=load('backup_full')

@contextmanager
def fixture():
    # Keep small test artifacts in the workspace; do not clean arbitrary temp paths.
    base=Path.cwd()/'tmp';base.mkdir(exist_ok=True)
    yield tempfile.mkdtemp(prefix='renamed_backup_test_',dir=base)

class Checks(unittest.TestCase):
    def test_bim_boundaries_unsorted_chromosomes(self):
        with fixture() as d:
            p=Path(d)/'test.bim'
            p.write_text('2 x 0 4 A C\n1 x 0 20 A C\n1 x 0 10 A C\n1 x 0 11 A C\n1 x 0 21 A C\n')
            q=[dict(chromosome='1',lo=10,hi=20),dict(chromosome='1',lo=11,hi=11),dict(chromosome='2',lo=5,hi=6)]
            self.assertEqual(bim.count_windows(p,q),([3,1,0],5))
    def test_backup_independence_and_corruption(self):
        with fixture() as d:
            project=Path(d)/'susie_mix';project.mkdir(); (project/'.hidden').write_text('hidden')
            (project/'results').mkdir();(project/'results/x.rds').write_bytes(b'fakeRDS')
            dest=Path(d)/'full_backup'
            calls=[]
            def copy(args,check):
                self.assertEqual(args[:4],['rsync','-a','--partial','--info=progress2'])
                shutil.copytree(project,dest/'susie_mix',dirs_exist_ok=True,symlinks=True)
                calls.append(args)
                return SimpleNamespace(returncode=24 if len(calls)==1 else 0)
            before=backup.inventory(project)
            with patch.object(backup.shutil,'which',return_value='/usr/bin/rsync'),patch.object(backup.subprocess,'run',side_effect=copy):
                backup.backup(project,dest)
            self.assertEqual(before,backup.inventory(project))
            self.assertTrue((dest/'susie_mix/.hidden').is_file())
            self.assertTrue((dest/'backup_metadata/BACKUP_VERIFIED.txt').is_file())
            self.assertFalse(os.path.samefile(project/'results/x.rds',dest/'susie_mix/results/x.rds'))
            with patch.object(backup,'digest',side_effect=AssertionError('Unchanged cached file rehashed')):
                self.assertTrue(backup.verify(project,dest/'susie_mix',dest/'backup_metadata')['exact_match'])
            (dest/'susie_mix/results/x.rds').write_bytes(b'badRDS!')
            result=backup.verify(project,dest/'susie_mix',dest/'backup_metadata')
            self.assertFalse(result['exact_match']);self.assertTrue(result['problems'])
    def test_backup_live_changes_are_not_verified(self):
        with fixture() as d:
            src=Path(d)/'source';dst=Path(d)/'copy';meta=Path(d)/'meta'
            for p in [src,dst,meta]:p.mkdir()
            (src/'x').write_text('original');shutil.copy2(src/'x',dst/'x')
            real=backup.digest
            def mutate(p):
                value=real(p)
                if p==src/'x':(src/'new').write_text('new data')
                return value
            with patch.object(backup,'digest',side_effect=mutate):result=backup.verify(src,dst,meta)
            self.assertFalse(result['exact_match']);self.assertIn('new',result['source_changes_during_verification'])
    def test_reject_backup_inside_source(self):
        with fixture() as d:
            with self.assertRaises(ValueError):backup.backup(Path(d),Path(d)/'backup')

if __name__=='__main__':unittest.main()
