import contextlib,csv,importlib.util,io,json,sys,tempfile,unittest,zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest import mock
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
PATCH=ROOT/'output/slider_subset_patch_2026-10-30'
spec=importlib.util.spec_from_file_location('subset_patch',PATCH/'install_subset.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
CONFIG=json.loads((PATCH/'patch.json').read_text())

class SubsetTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=ROOT/'tmp');self.root=Path(self.temp.name).resolve()
        assert self.root.parent==(ROOT/'tmp').resolve()
        self.target=self.root/'bundle';self.project=self.root/'project';self.state=self.root/'state'
        with zipfile.ZipFile(ROOT/'output/slider_inplace_rcc_2026-10-29.zip') as z:
            for n in z.namelist():
                p=self.target/('/'.join(n.split('/')[1:]));p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(z.read(n))
        self.reference=module.rows(self.target/'renamed_gene_map.csv')
        self.report=self.target/'checks/mapping_checks.csv';self.report.parent.mkdir()
        observed=[]
        for r in self.reference:
            x=dict(r,expression_matches='1',annotation_matches='1',issue='',status='PASS')
            for a,b in [('expression_id','actual_expression_id'),('expression_symbol','actual_expression_symbol'),('annotation_symbol','actual_annotation_symbol'),('chromosome','actual_chromosome'),('start','actual_start'),('end','actual_end')]:x[b]=r[a]
            if r['gene'] in CONFIG['deferred_genes']:x.update(actual_end=str(int(r['end'])+1),issue=module.ISSUE,status='FAIL')
            observed.append(x)
        self.observed=observed;self.save_report()
    def tearDown(self):
        assert self.root.parent==(ROOT/'tmp').resolve();self.temp.cleanup()
    def save_report(self):
        with self.report.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(self.observed[0]));w.writeheader();w.writerows(self.observed)
    def validate(self):return module.validate_report(self.target/'renamed_gene_map.csv',self.report,CONFIG)
    def test_fixed_646_34_partition(self):
        result=self.validate();self.assertEqual(result['accepted_ids'],646);self.assertEqual(result['deferred_ids'],34)
        kept=module.rows(PATCH/'replacement/renamed_gene_map.csv');excluded=module.rows(PATCH/'replacement/excluded_candidates.csv')
        original={r['gene']:r for r in self.reference}
        self.assertTrue(all(r==original[r['gene']] for r in kept))
        self.assertEqual(set(r['gene'] for r in excluded),set(CONFIG['deferred_genes']))
        self.assertEqual(len(set(r['expression_symbol'] for r in kept)),645)
    def test_extra_failures_are_not_silently_dropped(self):
        x=next(r for r in self.observed if r['status']=='PASS');x.update(actual_end=str(int(x['end'])+1),issue=module.ISSUE,status='FAIL');self.save_report()
        with self.assertRaisesRegex(RuntimeError,'Selection differs'):self.validate()
    def test_identity_or_ambiguity_changes_block_patch(self):
        for field,value in [('expression_matches','0'),('annotation_matches','2'),('actual_expression_id','ENSG_OTHER.1'),('actual_chromosome','X')]:
            old=self.observed[0][field];self.observed[0][field]=value;self.save_report()
            with self.assertRaises(RuntimeError):self.validate()
            self.observed[0][field]=old
    def test_duplicate_report_is_rejected(self):
        self.observed[-1]=self.observed[0].copy();self.save_report()
        with self.assertRaisesRegex(RuntimeError,'exactly'):self.validate()
    def test_interrupted_install_resumes_and_is_idempotent(self):
        def interrupt(rel):raise RuntimeError('interrupt')
        with self.assertRaisesRegex(RuntimeError,'interrupt'):module.install_files(self.target,self.state,CONFIG,after_write=interrupt)
        module.install_files(self.target,self.state,CONFIG);module.install_files(self.target,self.state,CONFIG)
        module.verify_complete(self.target,CONFIG['new_hashes'])
    def test_unexpected_file_change_prevents_any_replacement(self):
        before=(self.target/'renamed_gene_map.csv').read_bytes();(self.target/'initialize.R').write_text('unexpected')
        with self.assertRaisesRegex(RuntimeError,'modified'):module.install_files(self.target,self.state,CONFIG)
        self.assertEqual((self.target/'renamed_gene_map.csv').read_bytes(),before)
    def setup_project(self):
        d=self.project/'results_slide_prior_em/iteration_005/recovery/parallel';d.mkdir(parents=True)
        (d/'COMPLETE.rds').write_bytes(b'fixture certificate')
        p=self.project/'results_slide_prior_em/iteration_005/results/KEEP.rds';p.parent.mkdir();p.write_bytes(b'unchanged fit')
        return p
    def run_fixture(self,queue):
        @contextlib.contextmanager
        def lock(*args):yield
        storage=SimpleNamespace(maintenance_lock=lock,queue=lambda:queue,check_sources=lambda p:None,sync_dir=lambda p:None)
        fcntl=SimpleNamespace(flock=lambda *a:None,LOCK_EX=1,LOCK_NB=2)
        with mock.patch.dict(sys.modules,{'storage':storage,'fcntl':fcntl}):
            module.run(SimpleNamespace(project=self.project,target=self.target,state=self.state))
    def test_full_installer_preserves_005_and_backs_up_report(self):
        fit=self.setup_project();report=self.report.read_bytes();self.run_fixture([])
        self.assertEqual(fit.read_bytes(),b'unchanged fit')
        self.assertEqual((self.state/'mapping_checks.before_subset.csv').read_bytes(),report)
        self.assertEqual(module.rows(self.state/'original/renamed_gene_map.csv'),self.reference)
        self.run_fixture([]);self.assertTrue((self.state/'COMPLETE.json').exists())
    def test_active_jobs_block_install(self):
        self.setup_project()
        with self.assertRaisesRegex(RuntimeError,'still active'):self.run_fixture([['1','slide_next_gate','PENDING']])
        module.verify_complete(self.target,CONFIG['old_hashes'])
    def test_initialized_or_partial_006_is_never_redefined(self):
        self.setup_project();(self.project/'results_slide_prior_em/iteration_006').mkdir()
        with self.assertRaisesRegex(RuntimeError,'frozen cohort'):self.run_fixture([])
        module.verify_complete(self.target,CONFIG['old_hashes'])

if __name__=='__main__':unittest.main()
