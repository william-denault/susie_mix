import csv,importlib.util,json,os,sys,tempfile,unittest
from pathlib import Path
from unittest import mock
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import retention,storage,launch

class RetirementTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name).resolve()
        self.branch=self.root/'branch';self.branch.mkdir();self.source=self.root/'source';self.source.mkdir()
        storage.atomic_json(self.branch/'identity.json',dict(source=str(self.source)))
        self.state=self.branch/'state';self.state.mkdir()
        self.d=self.branch/'results_slide_prior_em/iteration_005';(self.d/'results').mkdir(parents=True)
        self.keep=self.branch/'results_slide_prior_em/iteration_006/results';self.keep.mkdir(parents=True)
        (self.keep/'KEEP.rds').write_bytes(b'keep')
        (self.d/'results/A.rds').write_bytes(b'old A');(self.d/'results/B.rds').write_bytes(b'old B')
        (self.state/'snapshot.rds').write_bytes(b'metadata')
        with (self.state/'directories.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=['directory']);w.writeheader();w.writerow(dict(directory=str(self.d)))
    def tearDown(self): self.tmp.cleanup()
    def plan(self): return retention.plan_files(self.branch,self.state,7)
    def run_plan(self,**kw): return retention.execute_plan(self.branch,self.state,7,**kw)
    def test_interrupted_deletion_resumes_and_keeps_previous(self):
        self.plan()
        def fail(p): raise RuntimeError('interrupt')
        with self.assertRaisesRegex(RuntimeError,'interrupt'): self.run_plan(after_unlink=fail)
        self.run_plan();self.run_plan()
        self.assertEqual(list((self.d/'results').iterdir()),[])
        self.assertEqual((self.keep/'KEEP.rds').read_bytes(),b'keep')
    def test_replacement_or_modified_old_file_blocks_all_deletion(self):
        self.plan();(self.d/'results/B.rds').write_bytes(b'changed')
        with self.assertRaisesRegex(RuntimeError,'changed'): self.run_plan()
        self.assertTrue((self.d/'results/A.rds').exists())
    def test_unexpected_file_blocks_all_deletion(self):
        self.plan();(self.d/'results/unexpected.txt').write_text('keep')
        with self.assertRaisesRegex(RuntimeError,'Unplanned'): self.run_plan()
        self.assertTrue((self.d/'results/A.rds').exists())
    def test_metadata_change_blocks_deletion(self):
        self.plan();(self.state/'snapshot.rds').write_bytes(b'changed')
        with self.assertRaisesRegex(RuntimeError,'metadata changed'): self.run_plan()
        self.assertTrue((self.d/'results/A.rds').exists())
    def test_current_directory_is_never_accepted(self):
        p=self.plan();p['directories']=[str(self.keep.parent)];storage.atomic_json(self.state/'plan.json',p)
        with self.assertRaisesRegex(RuntimeError,'another run'): self.run_plan()
        self.assertTrue((self.keep/'KEEP.rds').exists())
    def test_redirected_result_directory_refused(self):
        self.plan()
        with mock.patch.object(Path,'is_symlink',lambda p:str(p)==str(self.d/'results')):
            with self.assertRaisesRegex(RuntimeError,'redirected'): self.run_plan()
    def test_prior_and_pooled_files_are_not_retired(self):
        (self.d/'priors.csv').write_bytes(b'prior');(self.d/'recovery/parallel').mkdir(parents=True)
        (self.d/'recovery/parallel/pooled.rds').write_bytes(b'counts')
        (self.d/'recovery/parallel/COMPLETE.rds').write_bytes(b'cert')
        p=self.plan();self.run_plan();retention.archive_certificates(p);retention.archive_certificates(p)
        self.assertEqual((self.d/'priors.csv').read_bytes(),b'prior')
        self.assertEqual((self.d/'recovery/parallel/pooled.rds').read_bytes(),b'counts')
        self.assertEqual((self.d/'recovery/parallel/COMPLETE.before_rolling_retirement.rds').read_bytes(),b'cert')
    def test_live_workers_or_original_jobs_block_cleanup(self):
        with mock.patch.dict(os.environ,{'SLURM_JOB_ID':'123'}):
            own=['123','slide_expanded_control','RUNNING']
            for name in ('slide_expanded_fit','slide_parallel_control','slide_storage_compact'):
                with mock.patch.object(storage,'queue',return_value=[own,['124',name,'PENDING']]):
                    with self.assertRaisesRegex(RuntimeError,'still active'): retention.check_queue()
            with mock.patch.object(storage,'queue',return_value=[own,['125','simulation','RUNNING']]): retention.check_queue()
    def test_gate_waits_for_current_jobs_and_refuses_duplicates(self):
        rows=[['1','slide_storage_resume','PENDING'],['2','slide_parallel_fit','RUNNING'],['3','simulation','RUNNING']]
        self.assertEqual(len(launch.classify(rows)),2)
        for name in ('slide_roll_gate','slide_roll_setup','slide_expanded_control'):
            with self.assertRaisesRegex(RuntimeError,'duplicate'): launch.classify([['4',name,'PENDING']])

if __name__=='__main__': unittest.main()
