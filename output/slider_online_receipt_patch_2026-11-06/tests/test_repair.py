import ast, contextlib, hashlib, importlib.util, os, shutil, sys, unittest
from pathlib import Path
from unittest.mock import patch
HERE=Path(__file__).resolve().parents[1];ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
spec=importlib.util.spec_from_file_location('receipt_repair',HERE/'repair.py')
repair=importlib.util.module_from_spec(spec);spec.loader.exec_module(repair)
support=repair.support
prior=ROOT/'recovery/slider_online_submit_patch_2026-11-05/tests/test_repair.py'
spec=importlib.util.spec_from_file_location('submission_fixture',prior)
fixture=importlib.util.module_from_spec(spec);spec.loader.exec_module(fixture)

class Tests(unittest.TestCase):
    def setUp(self):
        base=fixture.RepairTests();base.setUp()
        self.__dict__.update({k:v for k,v in base.__dict__.items() if not k.startswith('_')})
        self.bundle=self.root/'receipt_patch';self.bundle.mkdir()
        shutil.copytree(HERE/'runtime',self.bundle/'runtime',ignore=shutil.ignore_patterns('__pycache__'))
        fixture.manifest(self.bundle/'runtime')
        shutil.copyfile(HERE/'baseline_SHA256SUMS',self.bundle/'baseline_SHA256SUMS')
        shutil.copyfile(HERE/'check_commit.R',self.bundle/'check_commit.R');fixture.manifest(self.bundle)
        for p in (ROOT/'output/slider_online_submit_patch_2026-11-05/runtime').rglob('*'):
            if p.is_file() and '__pycache__' not in p.parts:
                dest=self.branch/'code'/p.relative_to(ROOT/'output/slider_online_submit_patch_2026-11-05/runtime')
                dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        self.cfg['bundle_sha256']=hashlib.sha256((HERE/'baseline_SHA256SUMS').read_bytes()).hexdigest()
        support.atomic(self.branch/'config.json',self.cfg)
        frozen=support.read(self.branch/'frozen_inputs.json');frozen['config.json']=support.digest(self.branch/'config.json')
        support.atomic(self.branch/'frozen_inputs.json',frozen)
        p=self.branch/'plan/batch_03.txt';p.write_text('G5\nG6\n')
        pending=support.read(self.branch/'pending_submission.json');pending['job']='49853903'
        support.atomic(self.td/'submission.json',pending);(self.branch/'pending_submission.json').unlink()
        cmd=['Rscript','--vanilla',str(self.branch/'code/engine.R'),'commit',str(self.branch/'code'),str(self.branch),str(self.d)]
        support.atomic(self.branch/'ERROR.json',dict(job='49856781',error="Command '%s' returned non-zero exit status 1."%cmd))
        (self.branch/'logs').mkdir(exist_ok=True)
        (self.branch/'logs/control_49856781.err').write_text('Error in FUN(X[[i]], ...) : Missing/changed validated summary: T\n')
        self.state=self.branch/'repairs'/repair.TAG
        support.HERE=HERE;self.mb=support.runtime()
    def tearDown(self):fixture.RepairTests.tearDown(self)
    @contextlib.contextmanager
    def mocks(self):
        with patch.object(repair,'HERE',self.bundle),patch.object(support,'runtime',return_value=self.mb), \
             patch.object(self.mb,'lock',side_effect=lambda p:contextlib.nullcontext()),patch.object(self.mb,'queue',return_value=[]), \
             patch.object(repair.subprocess,'check_call') as preflight,patch.object(support,'resume') as resume:
            yield preflight,resume
    def test_metadata_and_all_results_preserved(self):
        before=support.protected_state(self.branch,self.mb)
        with self.mocks() as (preflight,resume):
            repair.repair(self.project,(2,2,2));preflight.assert_called_once();resume.assert_called_once()
        self.assertEqual(before,support.protected_state(self.branch,self.mb))
        self.assertEqual(support.read(self.td/'submission.json')['job'],'49853903')
        self.assertFalse((self.branch/'ERROR.json').exists());self.assertTrue((self.state/'ERROR.before_repair.json').exists())
        self.mb.verify_code(self.branch/'code');self.mb.verify_inputs(self.branch)
    def test_resumed_controller_commits_without_submitting_gene_jobs(self):
        with self.mocks():repair.repair(self.project,(2,2,2))
        with patch.dict(os.environ,SLURM_JOB_ID='999'),patch.object(self.mb,'lock',side_effect=lambda p:contextlib.nullcontext()), \
             patch.object(self.mb,'queue',return_value=[]),patch.object(self.mb,'submission_queue',return_value=[]), \
             patch.object(self.mb,'finish_source',return_value=True),patch.object(self.mb,'seed_run',return_value=True), \
             patch.object(self.mb,'rrun') as rrun,patch.object(self.mb,'submit',return_value='1000') as submit:
            self.mb.control(self.branch)
        self.assertEqual(submit.call_count,1);self.assertEqual(submit.call_args.args[1],'controller')
        self.assertIn('commit',[c.args[1] for c in rrun.call_args_list])
    def test_failed_preflight_leaves_error_and_code_unchanged(self):
        before=(self.branch/'code/SHA256SUMS').read_bytes()
        with self.mocks() as (preflight,resume):
            preflight.side_effect=RuntimeError('invalid summary')
            with self.assertRaisesRegex(RuntimeError,'invalid summary'):repair.repair(self.project,(2,2,2))
            resume.assert_not_called()
        self.assertTrue((self.branch/'ERROR.json').exists())
        self.assertEqual((self.branch/'code/SHA256SUMS').read_bytes(),before)
    def test_interrupted_swap_resumes_without_repeating_preflight(self):
        original=support.move
        def crash(src,dst,root):
            original(src,dst,root)
            if Path(src)==self.branch/'code':raise OSError('crash after code backup')
        with self.mocks(),patch.object(support,'move',side_effect=crash):
            with self.assertRaises(OSError):repair.repair(self.project,(2,2,2))
        with self.mocks() as (preflight,resume):
            repair.repair(self.project,(2,2,2));preflight.assert_not_called();resume.assert_called_once()
        self.mb.verify_code(self.branch/'code')
    def test_active_jobs_are_not_touched(self):
        with self.mocks(),patch.object(self.mb,'queue',return_value=[dict(name='slide_online_f_5af269')]):
            with self.assertRaisesRegex(RuntimeError,'jobs are active'):repair.repair(self.project,(2,2,2))
        self.assertTrue((self.branch/'ERROR.json').exists())
    def test_missing_receipt_blocks_repair(self):
        (self.d/'receipts/G6.json').unlink()
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'unvalidated or missing'):repair.repair(self.project,(2,2,2))
    def test_other_error_is_not_cleared(self):
        (self.branch/'logs/control_49856781.err').write_text('unrelated error')
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'literal-T error'):repair.repair(self.project,(2,2,2))
    def test_repeat_does_not_undo_later_progress(self):
        with self.mocks():repair.repair(self.project,(2,2,2))
        (self.branch/'progress.csv').write_bytes(b'stage\n3\n')
        with self.mocks() as (preflight,resume):
            repair.repair(self.project,(2,2,2));preflight.assert_not_called()
        self.assertEqual((self.branch/'progress.csv').read_bytes(),b'stage\n3\n')
    def test_python36(self):
        for p in (HERE/'repair.py',HERE/'repair_support.py',HERE/'build_bundle.py'):
            ast.parse(p.read_text(),feature_version=(3,6))

if __name__=='__main__':unittest.main()
