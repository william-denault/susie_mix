import ast, contextlib, hashlib, importlib.util, io, os, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
HERE=Path(__file__).resolve().parents[1];ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
import upgrade
s=upgrade.support

def manifest(p):
    files=sorted(f for f in p.rglob('*') if f.is_file() and f!=p/'SHA256SUMS' and '__pycache__' not in f.parts)
    (p/'SHA256SUMS').write_bytes(''.join(s.digest(f)+'  '+f.relative_to(p).as_posix()+'\n' for f in files).encode())

class UpgradeTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='transition_upgrade_',dir=str(ROOT/'tmp'))
        self.root=Path(self.temp.name).resolve();self.project=self.root/'project'
        self.branch=self.project/'workflows/slider_online_015_025';self.branch.mkdir(parents=True)
        self.bundle=self.root/'bundle';self.bundle.mkdir()
        shutil.copytree(HERE/'runtime',self.bundle/'runtime',ignore=shutil.ignore_patterns('__pycache__'));manifest(self.bundle/'runtime')
        shutil.copyfile(HERE/'baseline_SHA256SUMS',self.bundle/'baseline_SHA256SUMS');manifest(self.bundle)
        shutil.copytree(ROOT/'output/slider_online_queue_patch_2026-11-06/runtime',self.branch/'code',ignore=shutil.ignore_patterns('__pycache__'))
        cfg=dict(source=str(self.project),source_iteration=14,rho=.3,passes=6,batches=9,full_steps=5,
                 workers=298,id='abc123',bundle_sha256=s.digest(HERE/'baseline_SHA256SUMS'))
        s.atomic(self.branch/'config.json',cfg)
        for name in ('state.rds','model.rds','stages.csv','prior_history.csv','latest_priors.csv','original_prior_history.csv'):
            (self.branch/name).write_bytes(('protected '+name).encode())
        (self.branch/'progress.csv').write_text('stage\n7\n')
        (self.branch/'plan').mkdir();(self.branch/'plan/batch.txt').write_text('G\n')
        results=self.project/'results_slide_prior_em';results.mkdir()
        for name in ('prior_history.csv','prior_history_long.csv','online_latest_priors.csv'):(results/name).write_text('history')
        fit=results/'iteration_015/results/G.rds';fit.parent.mkdir(parents=True);fit.write_bytes(b'completed fit')
        self.fit=fit
        s.atomic(self.branch/'frozen_inputs.json',{n:s.digest(self.branch/n) for n in ('config.json','model.rds','stages.csv')})
        s.HERE=self.bundle;self.mb=s.runtime();self.state=upgrade.state_dir(self.branch);self.state.mkdir(parents=True)
        self.priority='100';self.workers=True;self.events=[]
        s.atomic(self.bundle/'launch.json',dict(gate='501',watchdog='502'))
    def tearDown(self):
        assert (ROOT/'tmp').resolve() in self.root.parents
        self.temp.cleanup()
    def queue(self):
        rows=[dict(job='101',array='101',name='slide_online_c_abc123',state='PENDING')]
        if self.workers:rows+=[dict(job='100_1',array='100',name='slide_online_f_abc123',state='RUNNING')]
        return rows
    def details(self,job):
        if job=='502':
            return dict(JobName='slide_cleanup_recover',UserId='test(1)',WorkDir=str(self.bundle),Command=str(self.bundle/'recover_job'),Dependency='afterany:501(unfulfilled)')
        self.assertEqual(job,'101')
        return dict(JobName='slide_online_c_abc123',UserId='test(1)',WorkDir=str(self.branch),
                    Command=str(self.branch/'code/controller_job'),JobState='PENDING',Priority=self.priority)
    def action(self,action,job):
        cmd=['scontrol',action,job]
        self.events.append(cmd)
        self.assertEqual(cmd[0],'scontrol')
        if cmd==['scontrol','hold','101']:self.priority='0'
        elif cmd==['scontrol','release','101']:self.priority='100'
        else:self.fail('Unexpected action: '+str(cmd))
    @contextlib.contextmanager
    def mocks(self):
        with patch.object(upgrade,'HERE',self.bundle),patch.object(upgrade,'runtime',return_value=self.mb), \
             patch.dict(os.environ,USER='test',SLURM_JOB_ID='501'), \
             patch.object(self.mb,'lock',side_effect=lambda p:contextlib.nullcontext()), \
             patch.object(self.mb,'queue',side_effect=self.queue),patch.object(upgrade,'details',side_effect=self.details), \
             patch.object(upgrade,'control_action',side_effect=self.action), \
             patch.object(upgrade.time,'sleep',side_effect=lambda n:setattr(self,'workers',False)):
            yield
    def held(self):
        self.priority='0';self.workers=False;s.atomic(self.state/'hold.json',dict(job='101',original_priority='100'))
    def test_live_gate_drains_preserves_and_releases(self):
        before=upgrade.protected(self.branch,self.mb);fit=s.digest(self.fit)
        with self.mocks():upgrade.gate(self.branch)
        self.assertEqual(self.events,[['scontrol','hold','101'],['scontrol','release','101']])
        self.assertEqual(before,upgrade.protected(self.branch,self.mb));self.assertEqual(fit,s.digest(self.fit))
        self.mb.verify_code(self.branch/'code');self.mb.verify_inputs(self.branch)
        self.assertTrue((self.state/'INSTALLED.json').exists())
    def test_no_install_while_workers_active(self):
        self.held();self.workers=True;before=(self.branch/'code/SHA256SUMS').read_bytes()
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'still active'):upgrade.install(self.branch,self.mb)
        self.assertEqual(before,(self.branch/'code/SHA256SUMS').read_bytes())
    def test_recovery_releases_hold_without_install_plan(self):
        self.held();self.workers=True;before=(self.branch/'code/SHA256SUMS').read_bytes()
        with self.mocks():upgrade.recover(self.branch,self.mb)
        self.assertEqual(self.priority,'100');self.assertEqual(before,(self.branch/'code/SHA256SUMS').read_bytes())
    def test_watchdog_finishes_interrupted_code_swap(self):
        self.held();original=s.move
        def interrupt(src,dst,root):
            original(src,dst,root)
            if Path(src)==self.branch/'code':raise OSError('node lost after code backup')
        with self.mocks(),patch.object(s,'move',side_effect=interrupt):
            with self.assertRaisesRegex(OSError,'node lost'):upgrade.install(self.branch,self.mb)
        self.assertFalse((self.branch/'code').exists())
        with self.mocks():upgrade.recover(self.branch,self.mb)
        self.mb.verify_code(self.branch/'code');self.mb.verify_inputs(self.branch);self.assertEqual(self.priority,'100')
    def test_watchdog_finishes_interrupted_configuration_update(self):
        self.held();original=upgrade.write
        def interrupt(path,data):
            original(path,data)
            if Path(path)==self.branch/'config.json':raise OSError('node lost after config write')
        with self.mocks(),patch.object(upgrade,'write',side_effect=interrupt):
            with self.assertRaisesRegex(OSError,'node lost'):upgrade.install(self.branch,self.mb)
        with self.mocks():upgrade.recover(self.branch,self.mb)
        self.mb.verify_code(self.branch/'code');self.mb.verify_inputs(self.branch);self.assertEqual(self.priority,'100')
    def test_install_failure_before_plan_restores_continuation(self):
        with self.mocks(),patch.object(upgrade,'protected',side_effect=RuntimeError('snapshot failed')):
            with self.assertRaisesRegex(RuntimeError,'snapshot failed'):upgrade.gate(self.branch)
        self.assertEqual(self.priority,'100');self.assertFalse((self.state/'PLAN.json').exists())
    def test_existing_error_is_never_cleared(self):
        s.atomic(self.branch/'ERROR.json',dict(error='another failure'))
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'Existing failure'):upgrade.gate(self.branch)
        self.assertFalse(self.events);self.assertTrue((self.branch/'ERROR.json').exists())
    def test_wrong_controller_identity_not_held(self):
        original=self.details
        def wrong(job):
            d=original(job)
            if job=='101':d['WorkDir']='/unrelated/project'
            return d
        with self.mocks(),patch.object(upgrade,'details',side_effect=wrong):
            with self.assertRaisesRegex(RuntimeError,'identity'):upgrade.gate(self.branch)
        self.assertFalse(self.events)
    def test_user_held_controller_not_released(self):
        self.priority='0'
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'already held'):upgrade.gate(self.branch)
        self.assertFalse(self.events);self.assertEqual(self.priority,'0')
    def test_no_confirmed_watchdog_no_hold(self):
        s.atomic(self.bundle/'launch.json',dict(gate='501'))
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'No confirmed recovery'):upgrade.gate(self.branch)
        self.assertFalse(self.events)
    def test_changed_analysis_blocks_install(self):
        self.held()
        with self.mocks():upgrade.install(self.branch,self.mb)
        (self.branch/'state.rds').write_text('changed checkpoint')
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'Analysis changed'):upgrade.install(self.branch,self.mb)
    def test_repeat_completed_upgrade_only_releases(self):
        with self.mocks():upgrade.gate(self.branch)
        before=s.digest(self.branch/'state.rds')
        with self.mocks(),patch.object(upgrade,'install') as install:
            upgrade.gate(self.branch);install.assert_not_called()
        self.assertEqual(before,s.digest(self.branch/'state.rds'))
    def test_launch_orders_held_gate_watchdog_then_release(self):
        (self.bundle/'launch.json').unlink();self.workers=False
        calls=[]
        def output(cmd,**kw):calls.append(cmd);return '501' if len(calls)==1 else '502'
        with self.mocks(),patch.object(upgrade.subprocess,'check_output',side_effect=output), \
             patch.object(upgrade,'control_action') as action:
            upgrade.launch()
        self.assertIn('--hold',calls[0]);self.assertIn('--dependency=afterany:501',calls[1])
        action.assert_called_once_with('release','501')
    def test_release_timeout_checks_state_before_another_action(self):
        self.held();calls=[]
        def accepted(action,job):
            self.action(action,job);calls.append((action,job))
            raise subprocess.TimeoutExpired(['scontrol',action,job],60)
        with self.mocks(),patch.object(upgrade,'control_action',side_effect=accepted):
            upgrade.release(self.branch,self.mb)
        self.assertEqual(calls,[('release','101')]);self.assertEqual(self.priority,'100')
    def test_release_transient_failure_retries_while_still_held(self):
        self.held();calls=[]
        def retry(action,job):
            calls.append((action,job))
            if len(calls)==1:raise subprocess.CalledProcessError(1,['scontrol',action,job],stderr='Socket timed out on send/recv operation')
            self.action(action,job)
        with self.mocks(),patch.object(upgrade,'control_action',side_effect=retry):
            upgrade.release(self.branch,self.mb)
        self.assertEqual(len(calls),2);self.assertEqual(self.priority,'100')
    def test_scontrol_read_timeout_is_retried(self):
        raw='JobName=fixture JobState=PENDING Priority=0'
        with patch.object(upgrade,'runtime',return_value=self.mb),patch.dict(os.environ,SLURM_JOB_ID='501'), \
             patch.object(self.mb.subprocess,'check_output',side_effect=[subprocess.TimeoutExpired(['scontrol'],60),raw]), \
             patch.object(self.mb.time,'sleep'),contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(upgrade.details('101')['Priority'],'0')
    def test_python36(self):
        for name in ('upgrade.py','support.py','build_bundle.py'):
            ast.parse((HERE/name).read_text(),feature_version=(3,6))

if __name__=='__main__':unittest.main()
