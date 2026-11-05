import ast, contextlib, importlib.util, os, tempfile, unittest
from pathlib import Path
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('mb',Path(__file__).resolve().parents[1]/'minibatch.py')
mb=importlib.util.module_from_spec(spec);spec.loader.exec_module(mb)

class SchedulerTests(unittest.TestCase):
    def setUp(self):
        scratch=Path('tmp').resolve();scratch.mkdir(exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(prefix='mb_scheduler_',dir=str(scratch));self.root=Path(self.temp.name).resolve()
        self.branch=self.root/'branch';self.branch.mkdir()
        self.cfg=dict(id='123456',workers=148,seed=20261104)
        mb.atomic(self.branch/'config.json',self.cfg)
    def tearDown(self):
        assert Path('tmp').resolve() in self.root.parents
        self.temp.cleanup()
    def test_plan(self):
        genes=['G%05d'%i for i in range(18263)]
        plan=mb.make_plan(self.branch,genes)
        self.assertEqual(len(plan),59)
        for k in range(6):
            batches=[(self.branch/r['genes_file']).read_text().splitlines() for r in plan[k*9:k*9+9]]
            self.assertEqual(sorted(map(len,batches)),[2029]*7+[2030]*2)
            self.assertEqual(sorted(sum(batches,[])),genes)
            self.assertEqual([r['rho'] for r in plan[k*9:k*9+9]],[.3]*9)
        self.assertNotEqual((self.branch/plan[0]['genes_file']).read_text(),(self.branch/plan[9]['genes_file']).read_text())
        for r in plan[54:]: self.assertEqual((self.branch/r['genes_file']).read_text().splitlines(),genes)
    def test_safe_windows(self):
        r=dict(name='slide_inplace_control',state='RUNNING',reason='node',start='Unknown')
        self.assertFalse(mb.safe_window([r],self.cfg,0))
        r.update(state='PENDING',reason='Priority',start='1970-01-01T00:10:00')
        self.assertFalse(mb.safe_window([r],self.cfg,0))
        r['reason']='BeginTime';self.assertTrue(mb.safe_window([r],self.cfg,0))
        self.assertFalse(mb.safe_window([r],self.cfg,550))
        self.assertTrue(mb.safe_window([dict(r,name=mb.names(self.cfg)[0])],self.cfg,550))
    def test_pins_and_retirement(self):
        src=self.root/'source.rds';src.write_bytes(b'immutable source')
        st=src.stat();expected=[dict(gene='A',source=str(src),size=st.st_size,mtime=st.st_mtime)]
        dest=self.branch/'seed/snapshot_001/results';mb.pin(expected,dest)
        self.assertTrue(os.path.samefile(src,dest/'A.rds'))
        src.unlink();self.assertEqual((dest/'A.rds').read_bytes(),b'immutable source')
        original=self.root/'original/results';original.mkdir(parents=True);(original/'B.rds').write_bytes(b'keep')
        with self.assertRaises(RuntimeError): mb.release_results(self.branch,original)
        plan=mb.make_plan(self.branch,['G%d'%i for i in range(9)])
        for r in plan[:27]:
            d=self.branch/r['directory'];(d/'results').mkdir(parents=True)
            (d/'results/G.rds').write_bytes(b'fit');(d/'history.csv').write_text('preserve')
        mb.retire(self.branch,plan[18],plan)
        self.assertFalse(list(dest.iterdir()))
        for r in plan[:9]: self.assertFalse(list((self.branch/r['directory']/'results').iterdir()))
        for r in plan[9:27]: self.assertTrue((self.branch/r['directory']/'results/G.rds').exists())
        for r in plan[:27]: self.assertEqual((self.branch/r['directory']/'history.csv').read_text(),'preserve')
        self.assertEqual((original/'B.rds').read_bytes(),b'keep')
    def test_budget_and_ambiguous_submission(self):
        rows=[dict(name='unrelated',job=str(i)) for i in range(299)]
        task=dict(directory=str(self.branch/'tasks/task_000001'),genes=['A','B'],stage=1)
        with self.assertRaises(RuntimeError): mb.submit(self.branch,'workers',rows,task)
        with patch.object(mb.subprocess,'check_output',return_value='12345\n') as call:
            mb.submit(self.branch,'controller',rows,dependency='12344')
            args=call.call_args[0][0]
            self.assertIn('--time=23:00:00',args);self.assertIn('--dependency=afterany:12344',args)
        with patch.object(mb.subprocess,'check_output',side_effect=mb.subprocess.TimeoutExpired('sbatch',10)):
            with self.assertRaises(mb.subprocess.TimeoutExpired):mb.submit(self.branch,'controller',[])
        self.assertTrue((self.branch/'pending_submission.json').exists())
        with self.assertRaises(RuntimeError):mb.submit(self.branch,'controller',[])
    def test_expanded_queue_and_assignments(self):
        lines='10_1|10|slide_mb_f_123456|RUNNING|node|Unknown\n10_2|10|slide_mb_f_123456|PENDING|Resources|N/A\n11|11|other|COMPLETED|node|N/A\n'
        with patch.dict(os.environ,USER='test'),patch.object(mb.subprocess,'check_output',return_value=lines) as call:
            rows=mb.queue();self.assertEqual(len(rows),2);self.assertIn('--array',call.call_args[0][0])
        p=self.branch/'tasks/task_000001/submission.json'
        mb.atomic(p,dict(job='10',task=dict(stage=1,genes=['A','B'])))
        active,attempts,n=mb.active_assignments(self.branch,rows,1)
        self.assertEqual(active,{'A','B'});self.assertEqual(attempts,dict(A=1,B=1));self.assertEqual(n,2)
        with self.assertRaises(RuntimeError): mb.active_assignments(self.branch,rows,2)
    def test_receipt_and_frozen_inputs(self):
        d=self.branch/'batch';d.mkdir();(d/'input.rds').write_bytes(b'prior')
        f=d/'A.rds';f.write_bytes(b'fit')
        mb.atomic(d/'receipts/A.json',dict(gene='A',input_sha256=mb.sha(d/'input.rds'),files={'A.rds':dict(stamp=mb.stamp(f))}))
        self.assertTrue(mb.is_done(self.branch,d,'A'));self.assertFalse(mb.is_done(self.branch,d,'B'))
        f.write_bytes(b'changed fit')
        with self.assertRaises(RuntimeError):mb.is_done(self.branch,d,'A')
        mb.atomic(self.branch/'frozen_inputs.json',{'config.json':mb.sha(self.branch/'config.json')})
        mb.verify_inputs(self.branch);mb.atomic(self.branch/'config.json',dict(changed=True))
        with self.assertRaises(RuntimeError):mb.verify_inputs(self.branch)
    def test_python36_syntax(self):
        ast.parse(Path(mb.__file__).read_text(),feature_version=(3,6))
    def test_successor_preserves_first_failure(self):
        original=dict(error='specific original validation failure',job='900')
        mb.atomic(self.branch/'ERROR.json',original)
        with patch.object(mb.sys,'argv',['minibatch.py','control',str(self.branch)]), \
             patch.object(mb,'control',side_effect=RuntimeError('stopped branch')):
            with self.assertRaises(RuntimeError):mb.main()
        self.assertEqual(mb.read(self.branch/'ERROR.json'),original)
    def test_original_controller_can_always_submit_its_successor(self):
        rows=[dict(name='other_worker') for _ in range(298)]
        rows.append(dict(name='slide_inplace_control',state='PENDING',reason='BeginTime',start='2099-01-01T00:00:00'))
        self.assertEqual(mb.submission_budget(rows,self.cfg),299)
        self.assertTrue(mb.safe_window(rows,self.cfg))
        with self.assertRaises(RuntimeError):mb.submit(self.branch,'controller',rows)
        rows.pop(0)
        with patch.object(mb.subprocess,'check_output',return_value='12345\n'):
            mb.submit(self.branch,'controller',rows)
        # At most 299 jobs after our submission: the original controller can
        # still add its own successor, then exit and release its current slot.
        self.assertLessEqual(len(rows)+1+1,300)
    def test_capture_race_defers_but_scientific_failure_stops(self):
        cfg=dict(self.cfg,source=str(self.root/'source'));mb.atomic(self.branch/'config.json',cfg)
        with patch.object(mb,'verify_source'),patch.object(mb,'rrun',side_effect=mb.subprocess.CalledProcessError(75,'Rscript')):
            self.assertIs(mb.seed_run(self.branch),False)
        self.assertFalse((self.branch/'seed/READY.json').exists())
        self.assertTrue((self.branch/'seed/snapshot_001/CAPTURE_FAILED.json').exists())
        with patch.object(mb,'verify_source'),patch.object(mb,'rrun',side_effect=mb.subprocess.CalledProcessError(1,'Rscript')):
            with self.assertRaises(mb.subprocess.CalledProcessError):mb.seed_run(self.branch)
        self.assertTrue((self.branch/'seed/snapshot_002/CAPTURE_FAILED.json').exists())
    def test_controller_reserves_successor_and_fills_only_free_slots(self):
        plan=mb.make_plan(self.branch,['G%05d'%i for i in range(1800)])
        mb.writecsv(self.branch/'progress.csv',[dict(stage=0)])
        d=self.branch/plan[0]['directory'];d.mkdir(parents=True);(d/'input.rds').write_bytes(b'prior')
        rows=[dict(job=str(i),name='other_worker',array=str(i)) for i in range(250)]
        rows.append(dict(job='900',name=mb.names(self.cfg)[0],array='900'))
        submissions=[]
        def submit(branch,kind,queue_rows,task=None,dependency=None):
            submissions.append((kind,task,dependency));n=len(task['genes']) if task else 1
            self.assertLessEqual(len(queue_rows)+n,300)
            if not task:rows.append(dict(job='901',name=mb.names(self.cfg)[0],array='901'))
            return '902' if task else '901'
        with patch.dict(os.environ,SLURM_JOB_ID='900'),patch.object(mb,'lock',lambda p:contextlib.nullcontext()), \
             patch.object(mb,'verify_code'),patch.object(mb,'seed_run'),patch.object(mb,'rrun'), \
             patch.object(mb,'queue',lambda:list(rows)),patch.object(mb,'submit',submit):
            mb.control(self.branch)
        self.assertEqual(submissions[0],('controller',None,'900'))
        self.assertEqual(submissions[1][0],'workers')
        self.assertEqual(len(submissions[1][1]['genes']),48)

if __name__=='__main__':unittest.main()
