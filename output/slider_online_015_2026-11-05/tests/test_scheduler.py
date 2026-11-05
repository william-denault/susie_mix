import ast, contextlib, sys, importlib.util, os, tempfile, unittest
from pathlib import Path
from unittest.mock import patch

base=Path(__file__).resolve().parents[1]
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
mb=module('online',base/'minibatch.py');handoff=module('handoff',base/'handoff.py')
class Tests(unittest.TestCase):
    def setUp(self):
        scratch=Path('tmp').resolve();scratch.mkdir(exist_ok=True)
        self.tmp=tempfile.TemporaryDirectory(prefix='online_scheduler_',dir=str(scratch));self.root=Path(self.tmp.name).resolve()
        self.project=self.root/'project';self.results=self.project/'results_slide_prior_em';self.results.mkdir(parents=True)
        self.branch=self.project/'workflows/slider_online_015_025';self.branch.mkdir(parents=True)
        self.cfg=dict(id='123456',workers=298,seed=20261104,source=str(self.project))
        mb.atomic(self.branch/'config.json',self.cfg);(self.branch/'logs').mkdir()
        self.plan=mb.make_plan(self.branch,['G%05d'%i for i in range(18263)])
    def tearDown(self):
        assert Path('tmp').resolve() in self.root.parents
        self.tmp.cleanup()
    def test_balanced_numbered_plan(self):
        self.assertEqual(len(self.plan),59)
        for p in range(6):
            rows=self.plan[p*9:p*9+9];self.assertEqual({r['iteration'] for r in rows},{15+p})
            batches=[(self.branch/r['genes_file']).read_text().splitlines() for r in rows]
            self.assertEqual(sorted(map(len,batches)),[2029]*7+[2030]*2)
            self.assertEqual(len(set(sum(batches,[]))),18263)
        self.assertEqual([r['iteration'] for r in self.plan[54:]],list(range(21,26)))
        self.assertNotEqual((self.branch/self.plan[0]['genes_file']).read_text(),(self.branch/self.plan[9]['genes_file']).read_text())
    def test_owned_directories_and_cleanup(self):
        (self.results/'iteration_014/results').mkdir(parents=True);(self.results/'iteration_014/results/G.rds').write_bytes(b'old fit')
        (self.results/'iteration_014/prior_history.csv').write_text('preserve')
        for i in (15,16,17):
            mb.assert_owned_iteration(self.branch,i);d=self.results/('iteration_%03d'%i)
            (d/'results').mkdir();(d/'results/G.rds').write_bytes(b'fit');(d/'prior_history.csv').write_text('preserve')
        pin=self.branch/'seed/snapshot_001/results';pin.mkdir(parents=True)
        src=self.results/'iteration_014/results/G.rds';st=src.stat()
        mb.pin([dict(gene='G',source=str(src),size=st.st_size,mtime=st.st_mtime)],pin)
        with patch.object(mb,'queue',return_value=[]):mb.retire(self.branch,self.plan[18],self.plan)
        self.assertFalse((pin/'G.rds').exists())
        for i in (14,15):self.assertFalse((self.results/('iteration_%03d'%i)/'results/G.rds').exists())
        for i in (16,17):self.assertTrue((self.results/('iteration_%03d'%i)/'results/G.rds').exists())
        for i in (14,15,16,17):self.assertEqual((self.results/('iteration_%03d'%i)/'prior_history.csv').read_text(),'preserve')
        (self.results/'iteration_018').mkdir()
        with self.assertRaises(RuntimeError):mb.assert_owned_iteration(self.branch,18)
    def test_no_cleanup_while_original_jobs_exist(self):
        with patch.object(mb,'queue',return_value=[dict(name='slide_inplace_fit')]):
            with self.assertRaises(RuntimeError):mb.retire(self.branch,self.plan[0],self.plan)
    def test_boundary_refuses_existing_015(self):
        d=self.results/'iteration_014';d.mkdir();(d/'settings.rds').write_bytes(b'fixture')
        mb.writecsv(self.results/'prior_history.csv',[dict(iteration=14)])
        mb.atomic(self.results/'rolling_cleanup_inplace/before_014/COMPLETE.json',dict(iteration=14))
        handoff.boundary(self.project,mb)
        (self.results/'iteration_015').mkdir()
        with self.assertRaises(RuntimeError):handoff.boundary(self.project,mb)
    def test_reserves_legacy_successor_slot(self):
        rows=[dict(name='other') for _ in range(298)]+[dict(name='slide_inplace_control',state='PENDING',reason='BeginTime',start='2099-01-01T00:00:00')]
        self.assertEqual(mb.submission_budget(rows,self.cfg),299)
        with self.assertRaises(RuntimeError):mb.submit(self.branch,'controller',rows)
    def test_old_jobs_gate_seed_and_do_not_touch_history(self):
        with patch.object(mb,'queue',return_value=[dict(name='slide_inplace_fit')]),patch.object(mb,'rrun') as r:
            self.assertFalse(mb.seed_run(self.branch));r.assert_not_called()
    def test_python36(self):
        for name in ('minibatch.py','handoff.py'):ast.parse((base/name).read_text(),feature_version=(3,6))
if __name__=='__main__':unittest.main()
