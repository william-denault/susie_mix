import contextlib, hashlib, importlib.util, os, shutil, sys, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
HERE=Path(__file__).resolve().parents[1];ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('seed_repair',HERE/'repair.py')
repair=importlib.util.module_from_spec(spec);spec.loader.exec_module(repair)

def manifest(root):
    files=sorted(p for p in root.rglob('*') if p.is_file() and p!=root/'SHA256SUMS' and '__pycache__' not in p.parts)
    (root/'SHA256SUMS').write_bytes(''.join(repair.digest(p)+'  '+p.relative_to(root).as_posix()+'\n' for p in files).encode())

class RepairTests(unittest.TestCase):
    def setUp(self):
        scratch=ROOT/'tmp';scratch.mkdir(exist_ok=True)
        self.tmp=tempfile.TemporaryDirectory(prefix='seed_patch_test_',dir=str(scratch));self.root=Path(self.tmp.name).resolve()
        self.bundle=self.root/'patch';self.bundle.mkdir();self.project=self.root/'project'
        self.branch=self.project/'workflows/slider_online_015_025';self.branch.mkdir(parents=True)
        self.results=self.project/'results_slide_prior_em';self.source=self.results/'iteration_014/results';self.source.mkdir(parents=True)
        self.first=self.results/'iteration_015/batches/batch_01';self.first.mkdir(parents=True)
        for sub in ('success','receipts','summaries'):(self.first/sub).mkdir()
        (self.results/'iteration_015/results').mkdir()
        self.seed=self.branch/'seed/snapshot_001';(self.seed/'results').mkdir(parents=True)
        shutil.copytree(HERE/'runtime',self.bundle/'runtime',ignore=shutil.ignore_patterns('__pycache__'))
        manifest(self.bundle/'runtime')
        old=ROOT/'output/slider_online_015_2026-11-05'
        shutil.copytree(old,self.branch/'code',ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copyfile(HERE/'baseline_SHA256SUMS',self.bundle/'baseline_SHA256SUMS')
        shutil.copyfile(HERE/'check_state.R',self.bundle/'check_state.R');manifest(self.bundle)
        self.cfg=dict(schema=3,source=str(self.project),branch=str(self.branch),source_iteration=14,
                      rho=.3,passes=6,batches=9,full_steps=5,workers=298,seed=20261104,id='abc123',poll_minutes=1,
                      bundle_sha256=hashlib.sha256((self.bundle/'baseline_SHA256SUMS').read_bytes()).hexdigest())
        repair.atomic(self.branch/'config.json',self.cfg)
        repair.atomic(self.branch/'ERROR.json',dict(error='Three worker attempts exhausted; no genes skipped: G1,G2',job='123'))
        repair.atomic(self.branch/'seed/READY.json',dict(snapshot=str(self.seed),source_iteration=14,genes=2))
        for rel in ('state.rds','model.rds','initial_priors.csv','original_prior_history.csv'):(self.branch/rel).write_bytes(b'unchanged '+rel.encode())
        (self.branch/'progress.csv').write_bytes(b'stage\n0\n')
        (self.branch/'stages.csv').write_bytes(b'stage,iteration\n1,15\n')
        (self.first/'input.rds').write_bytes(b'unchanged input')
        for name in ('prior_history.csv','prior_history_long.csv'):(self.results/name).write_bytes(b'unchanged history')
        rows=[]
        for gene in ('G1','G2'):
            p=self.source/(gene+'.rds');p.write_bytes(b'saved 014 '+gene.encode());st=p.stat()
            rows.append(dict(gene=gene,source=str(p),size=st.st_size,mtime=st.st_mtime,md5=hashlib.md5(p.read_bytes()).hexdigest()))
        mb=repair.module();mb.writecsv(self.seed/'expected.csv',rows)
        frozen={p.relative_to(self.branch).as_posix():repair.digest(p) for p in (self.branch/'config.json',self.branch/'state.rds',self.seed/'expected.csv')}
        repair.atomic(self.branch/'frozen_inputs.json',frozen)
        task=self.branch/'tasks/task_000001';task.mkdir(parents=True)
        repair.atomic(task/'task.json',dict(stage=1,genes=['G1','G2']))
        (task/'worker_123_1.err').write_text('Original failure log')
        self.state=self.branch/'repairs'/repair.TAG
        self.mb=mb
    def tearDown(self):
        assert (ROOT/'tmp').resolve() in self.root.parents
        self.tmp.cleanup()
    @contextlib.contextmanager
    def mocks(self):
        with patch.object(repair,'HERE',self.bundle),patch.object(repair,'module',return_value=self.mb), \
             patch.object(self.mb,'lock',side_effect=lambda p:contextlib.nullcontext()), \
             patch.object(self.mb,'queue',return_value=[]),patch.object(repair.subprocess,'check_call') as r, \
             patch.object(repair,'resume') as resume:
            yield r,resume
    def test_full_repair_preserves_state_and_archives_only_failed_tasks(self):
        oldstate=(self.branch/'state.rds').read_bytes();oldsource={p.name:repair.digest(p) for p in self.source.iterdir()}
        with self.mocks() as (r,resume):
            repair.repair(self.project,2)
            self.assertEqual(r.call_count,1);self.assertEqual(resume.call_count,1)
            repair.repair(self.project,2)
            self.assertEqual(r.call_count,1);self.assertEqual(resume.call_count,2)
        self.assertEqual(oldstate,(self.branch/'state.rds').read_bytes())
        self.assertEqual(oldsource,{p.name:repair.digest(p) for p in self.source.iterdir()})
        for p in self.source.iterdir():self.assertTrue(os.path.samefile(p,self.seed/'results'/p.name))
        self.assertFalse((self.branch/'ERROR.json').exists());self.assertEqual(list((self.branch/'tasks').iterdir()),[])
        self.assertEqual((self.state/'failed_tasks/task_000001/worker_123_1.err').read_text(),'Original failure log')
        self.mb.verify_code(self.branch/'code');self.mb.verify_inputs(self.branch)
    def test_interrupted_code_swap_can_resume(self):
        original=repair.move;fired=[False]
        def crash(source,dest,root):
            original(source,dest,root)
            if Path(source)==self.branch/'code' and not fired[0]:fired[0]=True;raise OSError('disconnect after backup')
        with self.mocks(),patch.object(repair,'move',side_effect=crash):
            with self.assertRaises(OSError):repair.repair(self.project,2)
        self.assertFalse((self.branch/'code').exists())
        with self.mocks():repair.repair(self.project,2)
        self.mb.verify_code(self.branch/'code');self.assertTrue((self.state/'COMPLETE.json').exists())
    def test_changed_source_stops_before_link_or_code_changes(self):
        (self.source/'G2.rds').write_bytes(b'changed')
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'Source result changed'):repair.repair(self.project,2)
        self.assertEqual(list((self.seed/'results').iterdir()),[]);self.assertTrue((self.branch/'ERROR.json').exists())
    def test_existing_results_or_committed_update_refuse_reset(self):
        (self.results/'iteration_015/results/G1.rds').write_bytes(b'result')
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'already contains'):repair.repair(self.project,2)
        self.assertTrue((self.branch/'tasks/task_000001').exists())
    def test_running_workers_block_repair(self):
        with self.mocks(),patch.object(self.mb,'queue',return_value=[dict(name='slide_online_f_abc123')]):
            with self.assertRaisesRegex(RuntimeError,'Slider jobs are active'):repair.repair(self.project,2)
    def test_R_preflight_failure_does_not_reset_attempts(self):
        with self.mocks() as (r,resume):
            r.side_effect=RuntimeError('wrong state')
            with self.assertRaisesRegex(RuntimeError,'wrong state'):repair.repair(self.project,2)
            resume.assert_not_called()
        self.assertTrue((self.branch/'tasks/task_000001').exists());self.assertEqual(list((self.seed/'results').iterdir()),[])
    def test_first_nine_batches_keep_seed_then_release_after_commit(self):
        mb=self.mb;mb.pin(repair.csvrows(self.seed/'expected.csv'),self.seed/'results')
        # The existing first online folder belongs to this synthetic workflow.
        mb.atomic(self.results/'iteration_015/ONLINE_OWNER.json',dict(workflow_id='abc123',iteration=15))
        for i in (16,17):mb.assert_owned_iteration(self.branch,i);(self.results/('iteration_%03d/results'%i)).mkdir()
        plan=[dict(stage=i+1,iteration=15) for i in range(9)]+[dict(stage=10,iteration=16)]
        with patch.object(mb,'queue',return_value=[]):
            for row in plan[:9]:
                mb.retire(self.branch,row,plan)
                self.assertEqual(len(list((self.seed/'results').glob('*.rds'))),2)
                self.assertEqual(len(list(self.source.glob('*.rds'))),2)
            (self.branch/'progress.csv').write_bytes(b'stage\n8\n')
            with self.assertRaisesRegex(RuntimeError,'entire first online pass'):mb.retire(self.branch,plan[9],plan)
            self.assertEqual(len(list((self.seed/'results').glob('*.rds'))),2)
            (self.branch/'progress.csv').write_bytes(b'stage\n9\n')
            mb.retire(self.branch,plan[9],plan)
            self.assertEqual(list((self.seed/'results').iterdir()),[])
            self.assertEqual(list(self.source.iterdir()),[])

if __name__=='__main__':unittest.main()
