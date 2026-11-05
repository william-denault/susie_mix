import ast, contextlib, hashlib, importlib.util, os, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
HERE=Path(__file__).resolve().parents[1];ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('submit_repair',HERE/'repair.py')
repair=importlib.util.module_from_spec(spec);spec.loader.exec_module(repair)
def manifest(root):
    files=sorted(p for p in root.rglob('*') if p.is_file() and p!=root/'SHA256SUMS' and '__pycache__' not in p.parts)
    (root/'SHA256SUMS').write_bytes(''.join(repair.digest(p)+'  '+p.relative_to(root).as_posix()+'\n' for p in files).encode())

class RepairTests(unittest.TestCase):
    def setUp(self):
        scratch=ROOT/'tmp';scratch.mkdir(exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(prefix='submit_repair_',dir=str(scratch));self.root=Path(self.temp.name).resolve()
        self.project=self.root/'project';self.branch=self.project/'workflows/slider_online_015_025';self.branch.mkdir(parents=True)
        self.results=self.project/'results_slide_prior_em';self.iteration=self.results/'iteration_015'
        (self.iteration/'results').mkdir(parents=True);(self.results/'iteration_014/results').mkdir(parents=True)
        (self.results/'iteration_014/results/ORIGINAL.rds').write_bytes(b'014 unchanged')
        self.bundle=self.root/'patch';self.bundle.mkdir()
        shutil.copytree(HERE/'runtime',self.bundle/'runtime',ignore=shutil.ignore_patterns('__pycache__'));manifest(self.bundle/'runtime')
        shutil.copyfile(HERE/'baseline_SHA256SUMS',self.bundle/'baseline_SHA256SUMS');manifest(self.bundle)
        shutil.copytree(ROOT/'output/slider_online_seed_patch_2026-11-05/runtime',self.branch/'code')
        self.mb=repair.runtime()
        self.cfg=dict(schema=3,source=str(self.project),branch=str(self.branch),source_iteration=14,id='5af269',
                      rho=.3,passes=6,batches=9,full_steps=5,workers=298,poll_minutes=1,seed=20261104,
                      bundle_sha256=hashlib.sha256((HERE/'baseline_SHA256SUMS').read_bytes()).hexdigest())
        repair.atomic(self.branch/'config.json',self.cfg)
        for name in ('state.rds','model.rds','prior_history.csv','initial_priors.csv','original_prior_history.csv','latest_priors.csv'):
            (self.branch/name).write_bytes(b'unchanged '+name.encode())
        for name in ('HANDOFF.json','seed/READY.json'):repair.atomic(self.branch/name,dict(kept=True))
        for name in ('prior_history.csv','prior_history_long.csv','online_latest_priors.csv'):
            (self.results/name).write_bytes(b'unchanged '+name.encode())
        (self.iteration/'prior_history.csv').write_bytes(b'unchanged history15')
        repair.atomic(self.iteration/'ONLINE_OWNER.json',dict(workflow_id='5af269',iteration=15))
        self.mb.writecsv(self.branch/'progress.csv',[dict(stage=2)])
        self.stages=[]
        for stage,genes in enumerate((['G1','G2'],['G3','G4'],['G5','G6','G7']),1):
            group='../../results_slide_prior_em/iteration_015';directory=group+'/batches/batch_%02d'%stage
            gf='plan/batch_%02d.txt'%stage;f=self.branch/gf;f.parent.mkdir(exist_ok=True);f.write_text('\n'.join(genes)+'\n')
            row=dict(stage=stage,phase='minibatch',pass_=1,batch=stage,iteration=15,directory=directory,genes_file=gf,group=group)
            row['pass']=row.pop('pass_');self.stages.append(row)
            d=self.branch/directory
            for sub in ('receipts','summaries','success'):(d/sub).mkdir(parents=True,exist_ok=True)
            (d/'input.rds').write_bytes(b'frozen input '+str(stage).encode())
            (d/'priors.csv').write_bytes(b'frozen priors '+str(stage).encode())
            if stage<=2:(d/'commit.rds').write_bytes(b'committed update '+str(stage).encode())
            for gene in genes:
                if gene=='G7':continue
                result=self.iteration/'results'/(gene+'.rds');result.write_bytes(b'fit '+gene.encode())
                (d/'summaries'/(gene+'.rds')).write_bytes(b'counts '+gene.encode())
                (d/'success'/(gene+'.csv')).write_bytes(b'success '+gene.encode())
                names=['../../results/'+gene+'.rds','summaries/'+gene+'.rds','success/'+gene+'.csv']
                self.mb.atomic(d/'receipts'/(gene+'.json'),dict(gene=gene,input_sha256=self.mb.sha(d/'input.rds'),
                    files={n:dict(stamp=self.mb.stamp(d/n)) for n in names}))
        self.mb.writecsv(self.branch/'stages.csv',self.stages)
        self.td=self.branch/'tasks/task_000038';self.td.mkdir(parents=True)
        self.d=self.branch/self.stages[2]['directory']
        self.task=dict(directory=str(self.td),genes=['G5','G6'],input_sha256=self.mb.sha(self.d/'input.rds'),stage=3,stage_dir=str(self.d))
        repair.atomic(self.td/'task.json',self.task)
        cmd=['sbatch','--parsable','--export=ALL','--time=23:00:00','--exclude=midway2-0003,midway2-0011',
             '--chdir='+str(self.branch),'--job-name=slide_online_f_5af269',
             '--output='+str(self.td/'worker_%A_%a.out'),'--error='+str(self.td/'worker_%A_%a.err'),
             '--array=1-2',str(self.branch/'code/worker_job'),str(self.branch),str(self.td/'task.json')]
        repair.atomic(self.branch/'pending_submission.json',dict(command=cmd,kind='workers',task=self.task,time=1234))
        repair.atomic(self.branch/'ERROR.json',dict(job='49853687',error=str(subprocess.TimeoutExpired(cmd,10))))
        for index in (1,2):(self.td/('worker_49853903_%d.err'%index)).write_bytes(b'')
        self.report=''.join('49853903_%d|slide_online_f_5af269|tester|COMPLETED|0:0|%s\n'%(i,self.branch) for i in (1,2))
        frozen={name:self.mb.sha(self.branch/name) for name in ('config.json','model.rds','stages.csv')}
        repair.atomic(self.branch/'frozen_inputs.json',frozen)
        self.state=self.branch/'repairs'/repair.TAG
    def tearDown(self):
        assert (ROOT/'tmp').resolve() in self.root.parents
        self.temp.cleanup()
    @contextlib.contextmanager
    def mocks(self):
        with patch.object(repair,'HERE',self.bundle),patch.object(repair,'runtime',return_value=self.mb), \
             patch.object(self.mb,'lock',side_effect=lambda p:contextlib.nullcontext()),patch.object(self.mb,'queue',return_value=[]), \
             patch.dict(os.environ,USER='tester'),patch.object(repair.subprocess,'check_output',return_value=self.report), \
             patch.object(repair,'resume') as resume:
            yield resume
    def test_reconcile_preserves_all_results_commits_and_task_attempts(self):
        before=repair.protected_state(self.branch,self.mb)
        with self.mocks() as resume:repair.repair(self.project,2);resume.assert_called_once()
        self.assertEqual(before,repair.protected_state(self.branch,self.mb))
        self.assertEqual((self.results/'iteration_014/results/ORIGINAL.rds').read_bytes(),b'014 unchanged')
        self.assertFalse((self.branch/'ERROR.json').exists());self.assertFalse((self.branch/'pending_submission.json').exists())
        self.assertEqual(repair.read(self.td/'submission.json')['job'],'49853903')
        self.assertTrue((self.td/'worker_49853903_1.err').exists())
        active,attempts,n=self.mb.active_assignments(self.branch,[],3)
        self.assertEqual(attempts,dict(G5=1,G6=1));self.assertEqual(n,0)
        missing=[g for g in ['G5','G6','G7'] if not self.mb.is_done(self.branch,self.d,g)]
        self.assertEqual(missing,['G7'])
        self.mb.verify_code(self.branch/'code');self.mb.verify_inputs(self.branch)
    def test_real_controller_excludes_completed_genes_and_commits(self):
        with self.mocks():repair.repair(self.project,2)
        submitted=[];r_calls=[]
        def submit(branch,kind,rows,task=None,dependency=None):submitted.append((kind,task));return '900'
        with patch.dict(os.environ,SLURM_JOB_ID='899'),patch.object(self.mb,'lock',side_effect=lambda p:contextlib.nullcontext()), \
             patch.object(self.mb,'queue',return_value=[]),patch.object(self.mb,'submission_queue',return_value=[]), \
             patch.object(self.mb,'finish_source',return_value=True),patch.object(self.mb,'seed_run',return_value=True), \
             patch.object(self.mb,'rrun',side_effect=lambda b,mode,*a:r_calls.append(mode)), \
             patch.object(self.mb,'submit',side_effect=submit):
            self.mb.control(self.branch)
        self.assertEqual([kind for kind,t in submitted],['controller','workers'])
        self.assertEqual(submitted[1][1]['genes'],['G7']);self.assertNotIn('commit',r_calls)
        self.assertEqual(int(self.mb.csvrows(self.branch/'progress.csv')[0]['stage']),2)
    def test_repeat_after_progress_does_not_roll_back(self):
        with self.mocks():repair.repair(self.project,2)
        (self.branch/'progress.csv').write_bytes(b'stage\n3\n')
        with self.mocks():repair.repair(self.project,2)
        self.assertEqual((self.branch/'progress.csv').read_bytes(),b'stage\n3\n')
    def test_interrupted_install_resumes(self):
        original=repair.move
        def crash(src,dst,root):
            original(src,dst,root)
            if Path(src)==self.branch/'code':raise OSError('crash after backup')
        with self.mocks(),patch.object(repair,'move',side_effect=crash):
            with self.assertRaises(OSError):repair.repair(self.project,2)
        self.assertFalse((self.branch/'code').exists())
        with self.mocks():repair.repair(self.project,2)
        self.assertTrue((self.state/'COMPLETE.json').exists())
    def test_interruption_after_pending_archive_resumes(self):
        original=repair.move
        def crash(src,dst,root):
            original(src,dst,root)
            if Path(src)==self.branch/'pending_submission.json':raise OSError('crash after archive')
        with self.mocks(),patch.object(repair,'move',side_effect=crash):
            with self.assertRaises(OSError):repair.repair(self.project,2)
        with self.mocks():repair.repair(self.project,2)
        self.assertFalse((self.branch/'ERROR.json').exists());self.assertTrue((self.state/'COMPLETE.json').exists())
    def test_active_worker_blocks_changes(self):
        with self.mocks(),patch.object(self.mb,'queue',return_value=[dict(name='slide_online_f_5af269')]):
            with self.assertRaisesRegex(RuntimeError,'jobs are active'):repair.repair(self.project,2)
        self.assertTrue((self.branch/'ERROR.json').exists());self.assertFalse((self.td/'submission.json').exists())
    def test_bad_accounting_blocks_changes(self):
        with self.mocks(),patch.object(repair.subprocess,'check_output',return_value=self.report.replace('COMPLETED','RUNNING')):
            with self.assertRaisesRegex(RuntimeError,'identity/state'):repair.repair(self.project,2)
        self.assertTrue((self.branch/'pending_submission.json').exists())
    def test_incomplete_accounting_blocks_changes(self):
        with self.mocks(),patch.object(repair.subprocess,'check_output',return_value=self.report.splitlines()[0]):
            with self.assertRaisesRegex(RuntimeError,'confirm all'):repair.repair(self.project,2)
    def test_changed_result_blocks_changes(self):
        (self.iteration/'results/G5.rds').write_bytes(b'changed output')
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'output changed'):repair.repair(self.project,2)
        self.assertFalse((self.state/'PLAN.json').exists())
    def test_different_failure_is_not_cleared(self):
        repair.atomic(self.branch/'ERROR.json',dict(error='different failure',job='49853687'))
        with self.mocks():
            with self.assertRaisesRegex(RuntimeError,'Different failure'):repair.repair(self.project,2)
    def test_submission_deadline_and_safe_window(self):
        with patch.object(self.mb.subprocess,'check_output',return_value='901\n') as call:
            (self.branch/'pending_submission.json').unlink()
            self.mb.submit(self.branch,'controller',[])
        self.assertEqual(call.call_args.kwargs['timeout'],120)
        import datetime
        now=1700000000
        row=dict(name='other_control',state='PENDING',reason='BeginTime',
                 start=datetime.datetime.fromtimestamp(now+180,datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%S'))
        self.assertFalse(self.mb.safe_window([row],self.cfg,now))
        row['start']=datetime.datetime.fromtimestamp(now+400,datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%S')
        self.assertTrue(self.mb.safe_window([row],self.cfg,now))
    def test_another_timeout_keeps_journal_and_never_blindly_resubmits(self):
        (self.branch/'pending_submission.json').unlink()
        with patch.object(self.mb.subprocess,'check_output',side_effect=subprocess.TimeoutExpired('sbatch',120)) as call:
            with self.assertRaises(subprocess.TimeoutExpired):self.mb.submit(self.branch,'workers',[],task=self.task)
            with self.assertRaisesRegex(RuntimeError,'Ambiguous previous'):self.mb.submit(self.branch,'workers',[],task=self.task)
        self.assertEqual(call.call_count,1);self.assertTrue((self.branch/'pending_submission.json').exists())
    def test_python36_syntax(self):
        for p in [HERE/'repair.py',HERE/'build_bundle.py',HERE/'runtime/minibatch.py']:
            ast.parse(p.read_text(),feature_version=(3,6))

if __name__=='__main__':unittest.main()
