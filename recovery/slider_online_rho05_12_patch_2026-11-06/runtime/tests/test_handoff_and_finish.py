import contextlib, importlib.util, os, sys, tempfile, unittest
from pathlib import Path
from unittest.mock import patch

BASE=Path(__file__).resolve().parents[1]
def module(name):
    spec=importlib.util.spec_from_file_location(name,BASE/(name+'.py'))
    obj=importlib.util.module_from_spec(spec);sys.modules[name]=obj;spec.loader.exec_module(obj);return obj
mb=module('minibatch');handoff=module('handoff')

class RecoveryTests(unittest.TestCase):
    def setUp(self):
        scratch=Path('tmp').resolve();scratch.mkdir(exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(prefix='online15_handoff_',dir=str(scratch))
        self.path=Path(self.temp.name).resolve();self.project=self.path/'project'
        self.root=self.project/'results_slide_prior_em';self.iteration=self.root/'iteration_014'
        self.iteration.mkdir(parents=True);(self.iteration/'settings.rds').write_bytes(b'settings')
        self.branch=self.project/'workflows/slider_online_015_025';self.branch.mkdir(parents=True)
        (self.branch/'logs').mkdir()
        self.cfg=dict(source=str(self.project),id='abc123',workers=298,seed=20261104)
        mb.atomic(self.branch/'config.json',self.cfg);mb.writecsv(self.root/'prior_history.csv',[dict(iteration=14)])
        mb.atomic(self.root/'rolling_cleanup_inplace/before_014/COMPLETE.json',dict(iteration=14))
        (self.root/'parallel_continuation_job_id.txt').write_text('456')
        self.worker=dict(job='455_1',array='455',name='slide_inplace_fit',state='RUNNING',reason='node',start='N/A')
        self.controller=dict(job='456',array='456',name='slide_inplace_control',state='PENDING',reason='BeginTime',start='N/A')
        self.rows=[dict(self.controller),dict(self.worker)];self.held=False;self.events=[]
    def tearDown(self):
        assert Path('tmp').resolve() in self.path.parents
        self.temp.cleanup()
    def details(self,job='456',state='PENDING',priority=None):
        return 'JobId=%s JobName=slide_inplace_control UserId=test(123) JobState=%s Priority=%s Command=%s latest WorkDir=/tmp'%(
            job,state,('0' if self.held else '100') if priority is None else priority,
            self.project/'slider_inplace_runtime_2026-10-29/job/resume_slide_prior_em')
    def output(self,cmd,**kwargs):
        self.events.append(cmd)
        self.assertEqual(cmd[:4],['scontrol','show','job','-o'])
        return self.details(cmd[-1])
    def action(self,cmd,**kwargs):
        self.events.append(cmd)
        if cmd==['scontrol','hold','456']:self.held=True
        elif cmd==['scancel','--state=PENDING','--name=slide_inplace_control','--user=test','456']:
            self.assertTrue(self.held);self.rows=[r for r in self.rows if r['job']!='456']
        else:self.fail('Unexpected mutation: '+str(cmd))
    @contextlib.contextmanager
    def mocks(self):
        with patch.dict(os.environ,USER='test',SLURM_JOB_ID='901'),patch.object(mb,'verify_source'), \
             patch.object(mb,'queue',side_effect=lambda:list(self.rows)), \
             patch.object(handoff.subprocess,'check_output',side_effect=self.output), \
             patch.object(handoff.subprocess,'check_call',side_effect=self.action):yield
    def test_handoff_accepts_unknown_or_past_estimates_and_preserves_workers(self):
        for start,reason in [('N/A','(BeginTime)'),('Unknown','(QOSMaxJobsPerUserLimit)'),('2000-01-01T00:00:00','Priority')]:
            with self.subTest(start=start):
                self.rows=[dict(self.controller,start=start,reason=reason),dict(self.worker)];self.held=False
                with self.mocks():handoff.begin(self.branch,mb)
                self.assertEqual(self.rows,[self.worker])
                state=mb.read(self.branch/'HANDOFF.json')
                self.assertEqual((state['source_iteration'],state['first_online_iteration'],state['last_full_iteration']),(14,15,25))
        holds=[x for x in self.events if x[:2]==['scontrol','hold']]
        cancels=[x for x in self.events if x[0]=='scancel']
        self.assertEqual(len(holds),3);self.assertEqual(len(cancels),3)
        self.assertTrue(all(x[-1]=='456' for x in cancels))
    def test_real_slurm_reason_format_is_normalized(self):
        for reason in ('(BeginTime)','BeginTime'):
            raw='456|456|slide_inplace_control   |PENDING|'+reason+'|2099-01-01T00:00:00\n'
            with patch.dict(os.environ,USER='test'),patch.object(mb.subprocess,'check_output',return_value=raw) as run:
                rows=mb.queue()
                self.assertEqual(rows[0]['reason'],'BeginTime')
                self.assertTrue(mb.safe_window(rows,self.cfg,now=0))
                self.assertIn('--format=%i|%F|%100j|%T|%r|%S',run.call_args[0][0])
    def test_no_cancel_if_hold_does_not_take_effect(self):
        with self.mocks(),patch.object(handoff,'inspect_controller',return_value={
            'JobState':'PENDING','Priority':'100'}):
            with self.assertRaisesRegex(RuntimeError,'not held'):handoff.begin(self.branch,mb)
        self.assertFalse(any(x[0]=='scancel' for x in self.events))
    def test_identity_failure_precedes_any_mutation(self):
        with self.mocks(),patch.object(handoff.subprocess,'check_output',return_value=self.details().replace('UserId=test','UserId=other')):
            with self.assertRaisesRegex(RuntimeError,'identity'):handoff.begin(self.branch,mb)
        self.assertFalse(any(x[:2]==['scontrol','hold'] or x[0]=='scancel' for x in self.events))
    def test_resume_after_cancel_does_not_cancel_again(self):
        real=mb.atomic
        def interrupted(path,obj):
            if Path(path).name=='HANDOFF.json':raise OSError('simulated disconnect')
            return real(path,obj)
        with self.mocks(),patch.object(mb,'atomic',side_effect=interrupted):
            with self.assertRaises(OSError):handoff.begin(self.branch,mb)
        with self.mocks():handoff.begin(self.branch,mb)
        self.assertEqual(len([x for x in self.events if x[0]=='scancel']),1)
    def test_running_race_is_released_not_cancelled(self):
        raced=[False]
        def output(cmd,**kwargs):
            self.events.append(cmd)
            return self.details(state='RUNNING' if raced[0] else 'PENDING')
        def action(cmd,**kwargs):
            self.events.append(cmd)
            if cmd[:2]==['scontrol','hold']:raced[0]=True
            elif cmd[:2]==['scontrol','release']:
                self.rows=[self.worker]
            else:self.fail('Running controller must not be cancelled')
        with self.mocks(),patch.object(handoff.subprocess,'check_output',side_effect=output), \
             patch.object(handoff.subprocess,'check_call',side_effect=action),patch.object(handoff.time,'sleep'):
            handoff.begin(self.branch,mb)
        self.assertFalse(any(x[0]=='scancel' for x in self.events))
    def test_existing_015_refused_before_hold(self):
        (self.root/'iteration_015').mkdir()
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'later iteration'):handoff.begin(self.branch,mb)
        self.assertEqual(self.events,[])
    def finish_case(self,status,workers=1,ambiguous=False):
        mb.atomic(self.branch/'HANDOFF.json',dict(source_iteration=14))
        saved=self.iteration/'results/ALREADY.rds';saved.parent.mkdir(exist_ok=True);saved.write_bytes(b'completed result')
        batch=self.iteration/'recovery/parallel/batches/batch_test';batch.mkdir(parents=True)
        if status=='complete':(self.iteration/'recovery/parallel/COMPLETE.rds').write_bytes(b'certified')
        rows=[dict(self.worker,job='455_'+str(i)) for i in range(workers)]+[
            dict(job='901',array='901',name=mb.names(self.cfg)[0],state='RUNNING',reason='node',start='N/A'),
            dict(job='902',array='902',name=mb.names(self.cfg)[0],state='PENDING',reason='Dependency',start='N/A')]
        calls=[]
        def execute(cmd,**kwargs):
            calls.append((cmd,kwargs))
            if cmd[0]=='sbatch':
                if ambiguous:raise mb.subprocess.TimeoutExpired(cmd,10)
                return '777\n'
            self.assertEqual(cmd[-5:-2],['plan',str(self.project),'14'])
            self.assertEqual(kwargs['env']['SUSIE_SLIDE_TARGET_ITERATION'],'14')
            if status=='batch':return 'batch\n14\n%s\n2\nfit\n'%batch
            return status+'\n14\n'
        with patch.dict(os.environ,USER='test',SLURM_JOB_ID='901'),patch.object(mb,'verify_source'), \
             patch.object(mb,'submission_queue',return_value=rows),patch.object(mb,'lock',side_effect=lambda p:contextlib.nullcontext()), \
             patch.object(mb.subprocess,'check_output',side_effect=execute):
            if ambiguous:
                with self.assertRaises(mb.subprocess.TimeoutExpired):mb.finish_source(self.branch)
            else:self.assertEqual(mb.finish_source(self.branch),status=='complete')
        self.assertEqual(saved.read_bytes(),b'completed result')
        return calls,batch
    def test_finish_reuses_saved_results_and_delegates_only_missing_tasks(self):
        calls,batch=self.finish_case('batch',workers=296)
        self.assertEqual(len(calls),2);cmd,kwargs=calls[-1]
        self.assertIn('--array=1-2',cmd);self.assertIn('--job-name=slide_inplace_fit',cmd)
        self.assertEqual(cmd[-2:],['14',str(batch)])
        self.assertEqual((batch/'array_job_id.txt').read_text().strip(),'777')
        self.assertFalse((self.branch/'pending_submission.json').exists())
    def test_full_budget_no_new_workers(self):
        calls,_=self.finish_case('waiting',workers=298)
        self.assertEqual(len(calls),1);self.assertEqual(calls[0][0][-2],'0')
    def test_certified_finish_unlocks_online_seed(self):
        calls,_=self.finish_case('complete',workers=0)
        self.assertEqual(len(calls),1);self.assertTrue((self.branch/'SOURCE_014_COMPLETE.json').exists())
    def test_ambiguous_submission_is_not_retried(self):
        self.finish_case('batch',workers=296,ambiguous=True)
        self.assertTrue((self.branch/'pending_submission.json').exists())
    def test_one_minute_continuation_has_dependency_and_budget(self):
        with patch.object(mb.subprocess,'check_output',return_value='903\n') as submit:
            mb.submit(self.branch,'controller',[],dependency='901')
            cmd=submit.call_args[0][0]
            self.assertIn('--begin=now+1minutes',cmd);self.assertIn('--dependency=afterany:901',cmd)
        rows=[dict(name='worker') for _ in range(300)]
        with self.assertRaises(RuntimeError):mb.submit(self.branch,'controller',rows)

if __name__=='__main__':unittest.main()
