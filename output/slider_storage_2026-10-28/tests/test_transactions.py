import csv,importlib.util,json,os,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.dont_write_bytecode=True
SOURCE=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('storage',SOURCE/'storage.py')
storage=importlib.util.module_from_spec(spec);spec.loader.exec_module(storage)

class Transactions(unittest.TestCase):
    def setUp(self):
        base=Path.cwd()/'tmp';base.mkdir(exist_ok=True)
        self.root=Path(tempfile.mkdtemp(prefix='slider_storage_test_',dir=str(base)))
        self.project=self.root/'project';self.project.mkdir()
        self.rel='results_slide_prior_em/iteration_004/results/GENE.rds'
        self.file=self.project/self.rel;self.file.parent.mkdir(parents=True)
        self.file.write_bytes(b'posterior:exact\n'+b'redundant call '*1000)
        self.old=self.file.read_bytes()
        self.ledger=storage.Ledger(self.root/'receipts.sqlite')
        self.stage=self.root/'scratch'
    def compact(self,source,dest,report):
        self.assertEqual(source.read_bytes(),self.old)
        dest.write_bytes(b'posterior:exact\n')
        report.write_text('removed_calls\n1\n')
    def run_one(self): return storage.transact(self.project,self.rel,self.stage,self.ledger,self.compact)
    def test_verified_replacement_and_repeat(self):
        r=self.run_one()
        self.assertEqual(self.file.read_bytes(),b'posterior:exact\n')
        self.assertGreater(r['old_bytes'],r['new_bytes'])
        self.assertEqual(self.run_one(),r)
        self.assertFalse((self.stage/'original.rds').exists())
    def test_interruption_after_original_vacated(self):
        real=storage.copy_verified
        def fail(source,dest,expected=None):
            if Path(dest).name.startswith('.storage_compact_'): raise OSError('Disk quota exceeded')
            return real(source,dest,expected)
        with patch.object(storage,'copy_verified',side_effect=fail):
            with self.assertRaisesRegex(OSError,'quota'): self.run_one()
        self.assertFalse(self.file.exists())
        self.assertEqual((self.stage/'original.rds').read_bytes(),self.old)
        self.run_one()
        self.assertEqual(self.file.read_bytes(),b'posterior:exact\n')
    def test_failed_verification_keeps_original(self):
        def bad(source,dest,report): raise RuntimeError('Scientific mismatch')
        with self.assertRaisesRegex(RuntimeError,'Scientific'):
            storage.transact(self.project,self.rel,self.stage,self.ledger,bad)
        self.assertEqual(self.file.read_bytes(),self.old)
    def test_unexpected_external_change_is_rejected(self):
        self.run_one();self.file.write_bytes(b'changed')
        with self.assertRaisesRegex(RuntimeError,'changed'): self.run_one()
    def test_corrupt_escrow_blocks_replacement(self):
        real=storage.copy_verified
        def fail(source,dest,expected=None):
            if Path(dest).name.startswith('.storage_compact_'): raise OSError('quota')
            return real(source,dest,expected)
        with patch.object(storage,'copy_verified',side_effect=fail):
            with self.assertRaises(OSError): self.run_one()
        (self.stage/'original.rds').write_bytes(b'corrupted')
        with self.assertRaisesRegex(RuntimeError,'Escrow'): self.run_one()
    def test_noop_does_not_rewrite(self):
        before=storage.stamp(self.file)
        def noop(source,dest,report):
            dest.write_bytes(b'not used');report.write_text('removed_calls\n0\n')
        r=storage.transact(self.project,self.rel,self.stage,self.ledger,noop)
        self.assertEqual(storage.stamp(self.file),before)
        self.assertEqual(r['old_md5'],r['new_md5'])
    def test_receipt_export(self):
        self.run_one();p=self.root/'receipts.csv';self.ledger.export(p)
        with p.open() as f: rows=list(csv.DictReader(f))
        self.assertEqual(rows[0]['relative'],self.rel)
        self.assertEqual(rows[0]['new_md5'],storage.hashes(self.file)['md5'])
    def test_interruption_after_replacement_before_commit(self):
        with patch.object(self.ledger,'put',side_effect=OSError('Interrupted receipt commit')):
            with self.assertRaisesRegex(OSError,'commit'): self.run_one()
        self.assertEqual(self.file.read_bytes(),b'posterior:exact\n')
        self.assertEqual((self.stage/'original.rds').read_bytes(),self.old)
        self.run_one()
        self.assertIsNotNone(self.ledger.get(self.rel))
    def test_hardlinked_source_is_not_vacated(self):
        alias=self.file.with_name('ALIAS.rds');os.link(self.file,alias)
        with self.assertRaisesRegex(RuntimeError,'Hard-linked'):self.run_one()
        self.assertEqual(self.file.read_bytes(),self.old)
        self.assertEqual(alias.read_bytes(),self.old)
    def test_inventory_does_not_touch_backup_or_later_directories(self):
        for name in ('iteration_000_backup','iteration_006'):
            p=self.project/'results_slide_prior_em'/name/'results'/'GENE.rds'
            p.parent.mkdir(parents=True);p.write_bytes(b'preserved')
        self.assertEqual(storage.inventory(self.project,5),[self.rel])

if __name__=='__main__':unittest.main()
