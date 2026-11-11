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
        shutil.copytree(ROOT/'output/slider_online_rho05_18_patch_2026-11-06/runtime',self.branch/'code',ignore=shutil.ignore_patterns('__pycache__'))
        cfg=dict(source=str(self.project),source_iteration=14,rho=.5,passes=6,batches=18,full_steps=5,
                 online_schedule=[dict(iteration=i,batches=9 if i<17 else 18,rho=.3 if i<17 else .5) for i in range(15,21)],
                 schedule_change=dict(tag='rho05_18_batches_2026-11-06',start_iteration=17),
                 workers=298,id='abc123',seed=20261104,bundle_sha256=s.digest(HERE/'baseline_SHA256SUMS'))
        s.atomic(self.branch/'config.json',cfg)
        for name in ('state.rds','model.rds','stages.csv','prior_history.csv','latest_priors.csv','original_prior_history.csv'):
            (self.branch/name).write_bytes(('protected '+name).encode())
        (self.branch/'progress.csv').write_text('stage\n92\n')
        self.genes=['G%05d'%i for i in range(18263)]
        s.HERE=self.bundle;self.mb=s.runtime();self.mb.make_plan(self.branch,self.genes)
        # Match the installed patch's nested future-plan namespace.
        plan=self.mb.csvrows(self.branch/'stages.csv')
        for row in plan:
            if row['phase']=='minibatch' and int(row['iteration'])>=17:
                src=self.branch/row['genes_file'];dst=self.branch/'plan/rho05_18'/src.name
                dst.parent.mkdir(exist_ok=True);shutil.copyfile(src,dst);row['genes_file']=dst.relative_to(self.branch).as_posix()
        self.mb.writecsv(self.branch/'stages.csv',plan)
        results=self.project/'results_slide_prior_em';results.mkdir()
        for name in ('prior_history.csv','prior_history_long.csv','online_latest_priors.csv'):(results/name).write_text('history')
        for iteration in range(15,24):self.mb.assert_owned_iteration(self.branch,iteration)
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
            return dict(JobName='slide_extend_recover',UserId='test(1)',WorkDir=str(self.bundle),Command=str(self.bundle/'recover_job'),Dependency='afterany:501(unfulfilled)')
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
             patch.object(upgrade,'inspect_checkpoint',side_effect=self.inspect), \
             patch.object(upgrade.time,'sleep',side_effect=lambda n:setattr(self,'workers',False)):
            yield
    def inspect(self,branch,state):
        (state/'checkpoint.csv').write_text('stage,genes\n92,18263\n')
        (state/'cohort.txt').write_text('\n'.join(self.genes)+'\n')
    def held(self):
        self.priority='0';self.workers=False;s.atomic(self.state/'hold.json',dict(job='101',original_priority='100'))
    def test_live_gate_drains_preserves_and_releases(self):
        before=upgrade.protected(self.branch,self.mb);fit=s.digest(self.fit)
        with self.mocks():upgrade.gate(self.branch)
        self.assertEqual(self.events,[['scontrol','hold','101'],['scontrol','release','101']])
        after=upgrade.protected(self.branch,self.mb);self.assertTrue(all(after[k]==v for k,v in before.items()));self.assertEqual(fit,s.digest(self.fit))
        self.mb.verify_code(self.branch/'code');self.mb.verify_inputs(self.branch)
        self.assertTrue((self.state/'INSTALLED.json').exists())
    def test_no_install_while_workers_active(self):
        self.held();self.workers=True;before=(self.branch/'code/SHA256SUMS').read_bytes()
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'still active'):upgrade.install(self.branch,self.mb)
        self.assertEqual(before,(self.branch/'code/SHA256SUMS').read_bytes())
    def test_previous_patch_watchdog_finishes_before_new_hold(self):
        waited=[False];original=self.queue
        def queue():
            rows=original()
            if not waited[0]:rows.append(dict(job='99',array='99',name='slide_rho12_recover',state='PENDING'))
            return rows
        def sleep(seconds):
            if not waited[0]:
                self.assertFalse(self.events);waited[0]=True
            self.workers=False
        with self.mocks(),patch.object(self.mb,'queue',side_effect=queue),patch.object(upgrade.time,'sleep',side_effect=sleep):
            upgrade.gate(self.branch)
        self.assertTrue(waited[0]);self.assertEqual(self.priority,'100')
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
            with self.assertRaisesRegex(RuntimeError,'Incomplete failure'):upgrade.gate(self.branch)
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
        for name in ('upgrade.py','support.py','replan.py','build_bundle.py','repair_submission.py'):
            ast.parse((HERE/name).read_text(),feature_version=(3,6))




    def test_watchdog_finishes_interrupted_stage_plan_update(self):
        self.held();original=upgrade.write
        def interrupt(path,data):
            original(path,data)
            if Path(path)==self.branch/'stages.csv':raise OSError('node lost after stages write')
        with self.mocks(),patch.object(upgrade,'write',side_effect=interrupt):
            with self.assertRaisesRegex(OSError,'node lost'):upgrade.install(self.branch,self.mb)
        with self.mocks():upgrade.recover(self.branch,self.mb)
        self.mb.verify_inputs(self.branch);self.assertEqual(len(self.mb.csvrows(self.branch/'stages.csv')),708)
        self.assertEqual(self.priority,'100')

    def test_unknown_owner_or_changed_cohort_refused(self):
        self.inspect(self.branch,self.state)
        (self.state/'cohort.txt').write_text('\n'.join(self.genes[:-1]+['UNKNOWN'])+'\n')
        with self.assertRaisesRegex(RuntimeError,'cohort differs'):upgrade.replan.prepare(self.branch,self.mb,self.state)
        self.inspect(self.branch,self.state)
        s.atomic(self.project/'results_slide_prior_em/iteration_016/ONLINE_OWNER.json',dict(workflow_id='foreign'))
        with self.assertRaisesRegex(RuntimeError,'ownership'):upgrade.replan.prepare(self.branch,self.mb,self.state)

    def test_resume_uses_upgraded_settings_and_never_rebuilds_plan(self):
        # Match a real launch identity, then install and resume from that code.
        import types
        cfg=s.read(self.branch/'config.json')
        cfg.update(schema=3,branch=str(self.branch),poll_minutes=1,
                   id=hashlib.sha256(str(self.branch).encode()).hexdigest()[:6])
        # Runtime job identity is mocked only for this scheduler-free resume test.
        s.atomic(self.branch/'config.json',cfg)
        frozen=s.read(self.branch/'frozen_inputs.json');frozen['config.json']=s.digest(self.branch/'config.json')
        s.atomic(self.branch/'frozen_inputs.json',frozen)
        for iteration in range(15,24):
            s.atomic(self.project/('results_slide_prior_em/iteration_%03d/ONLINE_OWNER.json'%iteration),dict(workflow_id=cfg['id']))
        self.inspect(self.branch,self.state)
        upgraded,stages,files,summary=upgrade.replan.prepare(self.branch,self.mb,self.state)
        upgraded['bundle_sha256']=s.digest(self.branch/'code/SHA256SUMS')
        s.atomic(self.branch/'config.json',upgraded);s.atomic(self.branch/'HANDOFF.json',dict(done=True))
        before=s.digest(self.branch/'stages.csv')
        args=types.SimpleNamespace(source=str(self.project),branch=str(self.branch),workers=298,resume=True,handoff=False)
        with self.mocks(),patch.object(self.mb,'__file__',str(self.branch/'code/minibatch.py')), \
             patch.object(self.mb,'verify_code'),patch.object(self.mb,'queue',return_value=[]), \
             patch.object(self.mb,'submission_queue',return_value=[]),patch.object(self.mb,'submit') as submit, \
             patch.object(self.mb,'make_plan') as make:
            self.mb.launch(args);submit.assert_called_once();make.assert_not_called()
        self.assertEqual(s.digest(self.branch/'stages.csv'),before)

    def test_new_plan_write_interruption_recovers_without_partial_install_release(self):
        self.held();original=upgrade.write;interrupted=[False]
        def interrupt(path,data):
            original(path,data)
            if Path(path).parent==self.branch/'plan/extend_024_038_2026-11-09' and not interrupted[0]:
                interrupted[0]=True;raise OSError('node lost after first future batch')
        with self.mocks(),patch.object(upgrade,'write',side_effect=interrupt):
            with self.assertRaisesRegex(OSError,'node lost'):upgrade.install(self.branch,self.mb)
        self.assertEqual(self.priority,'0');self.assertFalse((self.state/'INSTALLED.json').exists())
        with self.mocks():upgrade.recover(self.branch,self.mb)
        self.mb.verify_inputs(self.branch);self.assertEqual(self.priority,'100')

    def test_controller_continues_beyond_95_and_finishes_at_708(self):
        with self.mocks():upgrade.gate(self.branch)
        for progress in (95,707,708):
            (self.branch/'progress.csv').write_text('stage\n%d\n'%progress)
            row=self.mb.csvrows(self.branch/'stages.csv')[min(progress,707)]
            commands=[]
            def run(branch,mode,*args):
                commands.append(mode)
                if mode=='commit':(self.branch/'progress.csv').write_text('stage\n708\n')
            with self.mocks(),contextlib.ExitStack() as stack:
                for name,value in [('verify_code',None),('submission_queue',[]),('submit','999'),
                                   ('finish_source',True),('seed_run',True),('queue',[]),
                                   ('retire',None),('assert_owned_iteration',None),
                                   ('active_assignments',({}, {}, 0)),('is_done',True)]:
                    stack.enter_context(patch.object(self.mb,name,return_value=value))
                stack.enter_context(patch.object(self.mb,'rrun',side_effect=run))
                self.mb.control_once(self.branch)
            if progress==95:
                self.assertFalse((self.branch/'DONE.json').exists());self.assertIn('prepare',commands)
            else:
                self.assertEqual(s.read(self.branch/'DONE.json')['stages'],708)
                (self.branch/'DONE.json').unlink()

    def test_plan_61_batches_rhos_prefix_and_final_full_steps(self):
        import csv
        old=self.mb.csvrows(self.branch/'stages.csv');self.inspect(self.branch,self.state)
        cfg,data,files,summary=upgrade.replan.prepare(self.branch,self.mb,self.state)
        new=list(csv.DictReader(io.StringIO(data.decode())))
        self.assertEqual(new[:93],old[:93]);self.assertEqual(len(new),708)
        self.assertEqual(sorted(summary['batch_sizes']),[299]*37+[300]*24)
        for iteration in range(24,39):
            rows=[r for r in new if int(r['iteration'])==iteration]
            expected=.8 if iteration<=25 else .5 if iteration<=28 else .3 if iteration<=33 else 1.
            self.assertEqual(len(rows),61 if iteration<=33 else 1)
            self.assertEqual({float(r['rho']) for r in rows},{expected})
            selected=sum([files[r['genes_file']].decode().splitlines() for r in rows],[])
            self.assertEqual(len(selected),18263);self.assertEqual(set(selected),set(self.genes))
            if iteration<=33:
                original=self.mb.balanced(self.genes,61,20261104+iteration-14)
                self.assertEqual([set(files[r['genes_file']].decode().splitlines()) for r in rows],[set(b) for b in original])
        self.assertEqual((cfg['passes'],cfg['full_steps']),(16,8))
        self.assertEqual(summary['observed_runtimes'],0)

    def test_started_024_refused_even_without_results(self):
        self.mb.assert_owned_iteration(self.branch,24);self.inspect(self.branch,self.state)
        with self.assertRaisesRegex(RuntimeError,'already exists'):
            upgrade.replan.prepare(self.branch,self.mb,self.state)

    def test_submitted_024_task_refused(self):
        self.inspect(self.branch,self.state)
        s.atomic(self.branch/'tasks/task_000001/task.json',dict(stage=94))
        with self.assertRaisesRegex(RuntimeError,'task after 023'):
            upgrade.replan.prepare(self.branch,self.mb,self.state)

    def test_completed_023_is_supported(self):
        self.inspect(self.branch,self.state)
        (self.state/'checkpoint.csv').write_text('stage,genes\n93,18263\n')
        (self.branch/'progress.csv').write_text('stage\n93\n')
        _,_,_,summary=upgrade.replan.prepare(self.branch,self.mb,self.state)
        self.assertEqual(summary['checkpoint_stage'],93)

    def test_stopped_idle_upgrade_submits_once_and_watchdog_does_not_duplicate(self):
        self.workers=False
        with self.mocks(),patch.object(self.mb,'queue',return_value=[]), \
             patch.object(self.mb,'submission_queue',return_value=[]),patch.object(self.mb,'submit',return_value='777') as submit:
            upgrade.gate(self.branch);upgrade.recover(self.branch,self.mb)
        submit.assert_called_once();self.assertFalse(self.events)
        self.assertTrue((self.state/'INSTALLED.json').exists())
        self.assertEqual(s.read(self.state/'CONTINUED.json')['job'],'777')

    def test_other_failure_never_cleared(self):
        s.atomic(self.branch/'ERROR.json',dict(error='math problem',job='999'))
        s.atomic(self.branch/'pending_submission.json',dict(kind='workers',task={}))
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'Not the reported'):upgrade.gate(self.branch)
        self.assertTrue((self.branch/'ERROR.json').exists());self.assertFalse(self.events)

    def reported_submission(self):
        row=self.mb.csvrows(self.branch/'stages.csv')[92];d=self.branch/row['directory']
        d.mkdir(parents=True);(d/'input.rds').write_bytes(b'original frozen 023 input')
        td=self.branch/'tasks/task_001411';td.mkdir(parents=True)
        task=dict(directory=str(td),stage=93,stage_dir=str(d),genes=self.genes[:130],input_sha256=self.mb.sha(d/'input.rds'))
        s.atomic(td/'task.json',task)
        cmd=['sbatch','--parsable','--export=ALL','--time=23:00:00','--exclude=midway2-0003,midway2-0011',
             '--chdir='+str(self.branch),'--job-name=slide_online_f_abc123','--output='+str(td/'worker_%A_%a.out'),
             '--error='+str(td/'worker_%A_%a.err'),'--array=1-130',str(self.branch/'code/worker_job'),str(self.branch),str(td/'task.json')]
        s.atomic(self.branch/'pending_submission.json',dict(kind='workers',task=task,command=cmd))
        s.atomic(self.branch/'ERROR.json',dict(job='50013914',error=str(subprocess.CalledProcessError(1,cmd))))
        (self.branch/'logs').mkdir();(self.branch/'logs/control_50013914.err').write_text('sbatch: error: Batch job submission failed: Unexpected message received\n')
        outdir=self.project/'results_slide_prior_em/iteration_023/results';outdir.mkdir()
        for k,g in enumerate(task['genes'],1):
            (td/('worker_50014103_%d.err'%k)).write_text('')
            out=outdir/(g+'.rds');out.write_bytes(b'completed fit')
            s.atomic(d/'receipts'/(g+'.json'),dict(gene=g,input_sha256=task['input_sha256'],completed=200,
              files={str(out):dict(stamp=self.mb.stamp(out))}))
        return ''.join('50014103_%d|slide_online_f_abc123|test|COMPLETED|0:0|%s\n'%(k,self.branch) for k in range(1,131))

    def test_stopped_023_repair_and_schedule_upgrade_are_one_transaction(self):
        accounting=self.reported_submission();self.workers=False
        before=s.digest(self.branch/'state.rds')
        with self.mocks(),patch.object(self.mb,'queue',return_value=[]),patch.object(self.mb,'queue_output',return_value=accounting), \
             patch.object(self.mb,'submission_queue',return_value=[]),patch.object(self.mb,'submit',return_value='777') as submit:
            upgrade.gate(self.branch);upgrade.recover(self.branch,self.mb)
        submit.assert_called_once()
        self.assertFalse((self.branch/'ERROR.json').exists());self.assertFalse((self.branch/'pending_submission.json').exists())
        self.assertEqual(s.read(self.branch/'tasks/task_001411/submission.json')['job'],'50014103')
        self.assertEqual(s.digest(self.branch/'state.rds'),before)
        self.assertTrue((self.state/'ERROR.json.before_repair').exists())
        self.assertEqual(len(self.mb.csvrows(self.branch/'stages.csv')),708)

    def test_recovery_after_first_error_marker_is_archived(self):
        accounting=self.reported_submission();self.workers=False
        s.atomic(self.state/'idle.json',dict(time=0));original=upgrade.move
        def interrupt(src,dst,root):
            original(src,dst,root)
            if Path(src)==self.branch/'ERROR.json':raise OSError('lost node after archive')
        with self.mocks(),patch.object(self.mb,'queue',return_value=[]),patch.object(self.mb,'queue_output',return_value=accounting), \
             patch.object(upgrade,'move',side_effect=interrupt):
            with self.assertRaisesRegex(OSError,'lost node'):upgrade.install(self.branch,self.mb)
        self.assertFalse((self.state/'INSTALLED.json').exists())
        with self.mocks(),patch.object(self.mb,'queue',return_value=[]),patch.object(self.mb,'submission_queue',return_value=[]), \
             patch.object(self.mb,'submit',return_value='777'):
            upgrade.recover(self.branch,self.mb)
        self.assertFalse((self.branch/'pending_submission.json').exists());self.assertTrue((self.state/'INSTALLED.json').exists())

if __name__=='__main__':unittest.main()
