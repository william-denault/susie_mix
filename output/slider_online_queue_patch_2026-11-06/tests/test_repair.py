import ast, contextlib, hashlib, importlib.util, io, os, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path
from unittest.mock import patch

HERE=Path(__file__).resolve().parents[1];ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
import repair
support=repair.support

def manifest(path):
    paths=sorted(p for p in path.rglob('*') if p.is_file() and p!=path/'SHA256SUMS' and '__pycache__' not in p.parts)
    (path/'SHA256SUMS').write_bytes(''.join(support.digest(p)+'  '+p.relative_to(path).as_posix()+'\n' for p in paths).encode())

class RepairTests(unittest.TestCase):
    def setUp(self):
        scratch=ROOT/'tmp';scratch.mkdir(exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(prefix='queue_repair_',dir=str(scratch));self.root=Path(self.temp.name).resolve()
        self.project=self.root/'project';self.branch=self.project/'workflows/slider_online_015_025'
        self.branch.mkdir(parents=True);self.results=self.project/'results_slide_prior_em';self.results.mkdir()
        self.bundle=self.root/'bundle';self.bundle.mkdir()
        shutil.copytree(HERE/'runtime',self.bundle/'runtime',ignore=shutil.ignore_patterns('__pycache__'))
        manifest(self.bundle/'runtime')
        for name in ('baseline_SHA256SUMS','check_commit.R'):shutil.copyfile(HERE/name,self.bundle/name)
        manifest(self.bundle)
        shutil.copytree(ROOT/'output/slider_online_receipt_patch_2026-11-06/runtime',self.branch/'code',
                        ignore=shutil.ignore_patterns('__pycache__'))
        support.HERE=HERE;self.mb=support.runtime()
        self.cfg=dict(schema=3,source=str(self.project),branch=str(self.branch),source_iteration=14,
            rho=.3,passes=6,batches=9,full_steps=5,workers=298,poll_minutes=1,seed=20261104,id='abc123',
            bundle_sha256=hashlib.sha256((HERE/'baseline_SHA256SUMS').read_bytes()).hexdigest())
        self.mb.atomic(self.branch/'config.json',self.cfg)
        self.plan=self.mb.make_plan(self.branch,['G%02d'%i for i in range(18)])
        for name in ('state.rds','model.rds','prior_history.csv','initial_priors.csv','original_prior_history.csv','latest_priors.csv'):
            (self.branch/name).write_text('protected '+name)
        self.mb.writecsv(self.branch/'progress.csv',[dict(stage=5)])
        for name in ('HANDOFF.json','seed/READY.json'):self.mb.atomic(self.branch/name,dict(fixture=True))
        for name in ('prior_history.csv','prior_history_long.csv','online_latest_priors.csv'):
            (self.results/name).write_text('protected '+name)
        self.mb.assert_owned_iteration(self.branch,15)
        out=self.results/'iteration_015/results';out.mkdir()
        for i,row in enumerate(self.plan[:6],1):
            d=self.branch/row['directory']
            for sub in ('receipts','summaries','success'):(d/sub).mkdir(parents=True,exist_ok=True)
            (d/'input.rds').write_text('stage input '+str(i))
            if i<6:(d/'commit.rds').write_text('committed '+str(i))
            genes=(self.branch/row['genes_file']).read_text().splitlines()
            for gene in genes:
                result=(out/(gene+'.rds')).resolve();result.write_text('saved fit '+gene)
                summary=(d/'summaries'/(gene+'.rds'));summary.write_text('summary '+gene)
                success=(d/'success'/(gene+'.csv'));success.write_text('receipt '+gene)
                paths=[str(result),'summaries/'+gene+'.rds','success/'+gene+'.csv']
                self.mb.atomic(d/'receipts'/(gene+'.json'),dict(gene=gene,input_sha256=self.mb.sha(d/'input.rds'),
                    files={p:dict(stamp=self.mb.stamp(d/p)) for p in paths}))
        self.current=self.branch/self.plan[5]['directory']
        command=['squeue','--noheader','--array','--states=all','--user','wdenault','--format=%i|%F|%100j|%T|%r|%S']
        self.mb.atomic(self.branch/'ERROR.json',dict(job='49863377',error=str(subprocess.TimeoutExpired(command,15))))
        self.mb.atomic(self.branch/'frozen_inputs.json',{n:self.mb.sha(self.branch/n) for n in ('config.json','model.rds','stages.csv')})
        (self.branch/'logs').mkdir();self.state=self.branch/'repairs'/repair.TAG

    def tearDown(self):
        assert (ROOT/'tmp').resolve() in self.root.parents
        self.temp.cleanup()

    @contextlib.contextmanager
    def mocks(self):
        with patch.object(repair,'HERE',self.bundle),patch.object(support,'runtime',return_value=self.mb), \
             patch.object(self.mb,'lock',side_effect=lambda p:contextlib.nullcontext()),patch.object(self.mb,'queue',return_value=[]), \
             patch.object(repair.subprocess,'check_call') as preflight,patch.object(support,'resume') as resume:
            yield preflight,resume

    def test_saved_fits_five_updates_and_receipts_unchanged(self):
        before=support.protected_state(self.branch,self.mb)
        with self.mocks() as (preflight,resume):
            repair.repair(self.project,(2,)*6);preflight.assert_called_once();resume.assert_called_once()
        self.assertEqual(before,support.protected_state(self.branch,self.mb))
        self.assertFalse((self.branch/'ERROR.json').exists())
        self.assertTrue((self.state/'ERROR.before_repair.json').exists())
        self.mb.verify_code(self.branch/'code');self.mb.verify_inputs(self.branch)

    def test_controller_retries_at_each_read_without_replaying_commit_or_submission(self):
        with self.mocks():repair.repair(self.project,(2,)*6)
        calls=[];state={'attempt':False}
        def output(command,**kwargs):
            self.assertEqual(command[0],'squeue')
            if not state['attempt']:
                state['attempt']=True;raise subprocess.TimeoutExpired(command,60)
            state['attempt']=False;return ''
        with patch.dict(os.environ,USER='test',SLURM_JOB_ID='999'), \
             patch.object(self.mb,'lock',side_effect=lambda p:contextlib.nullcontext()), \
             patch.object(self.mb.subprocess,'check_output',side_effect=output) as query, \
             patch.object(self.mb.time,'sleep'),patch.object(self.mb,'finish_source',return_value=True), \
             patch.object(self.mb,'seed_run',return_value=True), \
             patch.object(self.mb,'rrun',side_effect=lambda b,m,*a:calls.append(m)), \
             patch.object(self.mb,'submit',return_value='1000') as submit,contextlib.redirect_stderr(io.StringIO()):
            self.mb.control(self.branch)
        self.assertEqual(calls.count('commit'),1)
        self.assertEqual(submit.call_count,1)
        self.assertEqual(submit.call_args.args[1],'controller')
        self.assertGreaterEqual(query.call_count,8)
        self.assertFalse((self.branch/'ERROR.json').exists())

    def test_worker_submission_waits_for_fresh_queue_after_timeout(self):
        with self.mocks():repair.repair(self.project,(2,)*6)
        # An unfinished next batch: only one fake gene needs a worker.
        gene=(self.branch/self.plan[5]['genes_file']).read_text().splitlines()[0]
        (self.current/'receipts'/(gene+'.json')).unlink()
        calls=[];pending=[False]
        def output(cmd,**kw):
            self.assertEqual(cmd[0],'squeue')
            if not pending[0]:pending[0]=True;raise subprocess.TimeoutExpired(cmd,60)
            pending[0]=False;return ''
        with patch.dict(os.environ,USER='test',SLURM_JOB_ID='999'), \
             patch.object(self.mb,'lock',side_effect=lambda p:contextlib.nullcontext()), \
             patch.object(self.mb.subprocess,'check_output',side_effect=output),patch.object(self.mb.time,'sleep'), \
             patch.object(self.mb,'finish_source',return_value=True),patch.object(self.mb,'seed_run',return_value=True), \
             patch.object(self.mb,'rrun',side_effect=lambda b,m,*a:calls.append(m)), \
             patch.object(self.mb,'submit',return_value='1000') as submit,contextlib.redirect_stderr(io.StringIO()):
            self.mb.control(self.branch)
        self.assertNotIn('commit',calls)
        self.assertEqual([c.args[1] for c in submit.call_args_list],['controller','workers'])
        self.assertEqual(submit.call_args.kwargs['task']['genes'],[gene])

    def test_failed_preflight_keeps_error_and_old_code(self):
        before=(self.branch/'code/SHA256SUMS').read_bytes()
        with self.mocks() as (preflight,resume):
            preflight.side_effect=RuntimeError('invalid summary')
            with self.assertRaisesRegex(RuntimeError,'invalid summary'):repair.repair(self.project,(2,)*6)
            resume.assert_not_called()
        self.assertTrue((self.branch/'ERROR.json').exists())
        self.assertEqual((self.branch/'code/SHA256SUMS').read_bytes(),before)

    def test_interrupted_install_is_resumable(self):
        original=support.move
        def crash(src,dst,root):
            original(src,dst,root)
            if Path(src)==self.branch/'code':raise OSError('interrupted swap')
        with self.mocks(),patch.object(support,'move',side_effect=crash):
            with self.assertRaisesRegex(OSError,'interrupted swap'):repair.repair(self.project,(2,)*6)
        with self.mocks() as (preflight,resume):
            repair.repair(self.project,(2,)*6);preflight.assert_not_called();resume.assert_called_once()
        self.mb.verify_code(self.branch/'code')

    def test_active_jobs_are_left_alone(self):
        with self.mocks(),patch.object(self.mb,'queue',return_value=[dict(name='slide_online_c_abc123')]):
            with self.assertRaisesRegex(RuntimeError,'jobs are active'):repair.repair(self.project,(2,)*6)
        self.assertTrue((self.branch/'ERROR.json').exists())

    def test_missing_receipt_blocks_repair(self):
        next((self.current/'receipts').glob('*.json')).unlink()
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'missing validated receipts'):repair.repair(self.project,(2,)*6)

    def test_changed_fit_blocks_repair(self):
        next((self.results/'iteration_015/results').glob('*.rds')).write_text('changed result')
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'Completed output changed'):repair.repair(self.project,(2,)*6)

    def test_other_error_is_not_cleared(self):
        self.mb.atomic(self.branch/'ERROR.json',dict(job='49863377',error='actual fitting error'))
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'Not the reported squeue timeout'):repair.repair(self.project,(2,)*6)

    def test_ambiguous_submission_is_not_retried(self):
        self.mb.atomic(self.branch/'pending_submission.json',dict(kind='workers'))
        with self.mocks() as (preflight,resume):
            with self.assertRaisesRegex(RuntimeError,'Unresolved submission'):repair.repair(self.project,(2,)*6)
            preflight.assert_not_called();resume.assert_not_called()

    def test_repeat_does_not_roll_back_progress(self):
        with self.mocks():repair.repair(self.project,(2,)*6)
        self.mb.writecsv(self.branch/'progress.csv',[dict(stage=6)])
        with self.mocks() as (preflight,resume):
            repair.repair(self.project,(2,)*6);preflight.assert_not_called()
        self.assertEqual(int(self.mb.csvrows(self.branch/'progress.csv')[0]['stage']),6)

    def test_python36(self):
        for p in (HERE/'repair.py',HERE/'repair_support.py',HERE/'build_bundle.py',HERE/'runtime/minibatch.py'):
            ast.parse(p.read_text(),feature_version=(3,6))

if __name__=='__main__':unittest.main()
