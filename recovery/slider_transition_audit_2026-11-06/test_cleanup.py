import hashlib, importlib.util, tempfile, unittest
from pathlib import Path
from unittest.mock import patch

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('cleanup_candidate',HERE/'candidate_runtime/minibatch.py')
mb=importlib.util.module_from_spec(spec);spec.loader.exec_module(mb)

class CleanupTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='transition_cleanup_',dir=str(ROOT/'tmp'))
        self.project=Path(self.temp.name).resolve();self.branch=self.project/'workflows/slider_online_015_025';self.branch.mkdir(parents=True)
        self.root=self.project/'results_slide_prior_em';self.root.mkdir()
        mb.atomic(self.branch/'config.json',dict(source=str(self.project),source_iteration=14,id='fixture',workers=298,seed=20261104))
        self.plan=mb.make_plan(self.branch,['G%05d'%i for i in range(18263)])
        mb.writecsv(self.branch/'progress.csv',[dict(stage=9)])
        self.dirs={}
        for i in (14,15,16):
            if i>=15:mb.assert_owned_iteration(self.branch,i)
            d=self.root/('iteration_%03d'%i);(d/'results').mkdir(parents=True)
            (d/'results/G.rds').write_bytes(b'saved gene fit')
            (d/'prior_history.csv').write_bytes(b'prior history')
            (d/'pooled_counts.csv').write_bytes(b'counts')
            self.dirs[i]=d
        self.seed=self.branch/'seed/snapshot_001/results';self.seed.mkdir(parents=True)
        (self.seed/'G.rds').write_bytes(b'initial warm start')
    def tearDown(self):
        assert (ROOT/'tmp').resolve() in self.project.parents
        self.temp.cleanup()
    def retire(self,index):
        with patch.object(mb,'queue',return_value=[]):mb.retire(self.branch,self.plan[index-1],self.plan)
    def test_full_18263_gene_plan_once_per_pass(self):
        for offset in range(0,54,9):
            groups=[(self.branch/r['genes_file']).read_text().splitlines() for r in self.plan[offset:offset+9]]
            self.assertEqual(list(map(len,groups)),[2030,2030]+[2029]*7)
            self.assertEqual(len(set(sum(groups,[]))),18263)
        self.assertEqual([r['iteration'] for r in self.plan[54:]],list(range(21,26)))
    def test_source_kept_until_ninth_batch_committed(self):
        mb.writecsv(self.branch/'progress.csv',[dict(stage=8)])
        with self.assertRaisesRegex(RuntimeError,'entire first online pass'):self.retire(10)
        self.assertTrue((self.seed/'G.rds').exists());self.assertTrue((self.dirs[14]/'results/G.rds').exists())
    def test_recognized_atomic_orphans_retire_only_with_old_results(self):
        for i in (14,15,16):(self.dirs[i]/'results/.em_write_1a2b3c').write_bytes(b'interrupted save')
        self.retire(10);self.retire(10)
        self.assertFalse(list((self.dirs[14]/'results').iterdir()))
        self.assertFalse(list(self.seed.iterdir()))
        for i in (15,16):
            self.assertTrue((self.dirs[i]/'results/G.rds').exists())
            self.assertTrue((self.dirs[i]/'results/.em_write_1a2b3c').exists())
        for i in (14,15,16):
            self.assertEqual((self.dirs[i]/'prior_history.csv').read_bytes(),b'prior history')
            self.assertEqual((self.dirs[i]/'pooled_counts.csv').read_bytes(),b'counts')
    def test_unknown_file_blocks_deletion_before_first_unlink(self):
        for name in ('notes.csv','.em_write_history.csv','unexpected.txt'):
            with self.subTest(name=name):
                f=self.dirs[14]/'results'/name;f.write_bytes(b'keep')
                with self.assertRaisesRegex(RuntimeError,'Unexpected object'):self.retire(10)
                self.assertTrue((self.dirs[14]/'results/G.rds').exists());f.unlink()
    def test_directory_blocks_deletion(self):
        (self.dirs[14]/'results/.em_write_ab12').mkdir()
        with self.assertRaisesRegex(RuntimeError,'Unexpected object'):self.retire(10)
        self.assertTrue((self.dirs[14]/'results/G.rds').exists())
    def test_cleanup_interrupted_then_resumed(self):
        d=self.dirs[14]/'results'
        for name in ('H.rds','.em_write_b123'): (d/name).write_bytes(b'old fit')
        actual=Path.unlink;count=[0]
        def stop(p,*a,**kw):
            actual(p,*a,**kw)
            if p.parent==d:
                count[0]+=1
                if count[0]==1:raise KeyboardInterrupt('simulated kill')
        with patch.object(Path,'unlink',stop):
            with self.assertRaises(KeyboardInterrupt):self.retire(10)
        self.retire(10)
        self.assertFalse(list(d.iterdir()))
        self.assertTrue((self.dirs[15]/'results/G.rds').exists())
        self.assertTrue((self.dirs[16]/'results/G.rds').exists())
    def test_legacy_jobs_block_all_cleanup(self):
        with patch.object(mb,'queue',return_value=[dict(name='slide_inplace_fit')]):
            with self.assertRaisesRegex(RuntimeError,'cleanup is blocked'):mb.retire(self.branch,self.plan[9],self.plan)
        self.assertTrue((self.seed/'G.rds').exists())
        self.assertTrue((self.dirs[14]/'results/G.rds').exists())
    def test_final_controller_exits_without_submitting(self):
        mb.atomic(self.branch/'DONE.json',dict(stages=59))
        with patch.object(mb,'verify_code'),patch.object(mb,'submit') as submit:
            mb.control(self.branch);submit.assert_not_called()

if __name__=='__main__':unittest.main()
