import contextlib, importlib.util, os, subprocess, sys, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
HERE=Path(__file__).resolve().parents[1];ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE));sys.path.insert(0,str(HERE/'runtime'))
import support, repair_submission, runtime_order, submission_receipts
support.HERE=HERE
mb=support.runtime()


class ReceiptsOrderTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='transition_receipt_',dir=str(ROOT/'tmp'))
        self.root=Path(self.temp.name).resolve();self.branch=self.root/'workflows/slider_online_015_025'
        self.branch.mkdir(parents=True);self.genes=['G%03d'%i for i in range(130)]
        mb.atomic(self.branch/'config.json',dict(id='abc123',source=str(self.root),workers=298))
        self.d=self.root/'results_slide_prior_em/iteration_023/batches/full_step'
        self.d.mkdir(parents=True);(self.d/'input.rds').write_bytes(b'fixed input')
        self.td=self.branch/'tasks/task_001411';self.td.mkdir(parents=True)
        self.task=dict(directory=str(self.td),stage=93,stage_dir=str(self.d),genes=self.genes,input_sha256=mb.sha(self.d/'input.rds'))
        mb.atomic(self.td/'task.json',self.task)
        command=['sbatch','--parsable','--export=ALL','--time=23:00:00','--exclude=midway2-0003,midway2-0011',
          '--chdir='+str(self.branch),'--job-name=slide_online_f_abc123','--output='+str(self.td/'worker_%A_%a.out'),
          '--error='+str(self.td/'worker_%A_%a.err'),'--array=1-130',str(self.branch/'code/worker_job'),str(self.branch),str(self.td/'task.json')]
        self.pending=dict(kind='workers',task=self.task,command=command,time=1794252795.721)
        mb.atomic(self.branch/'pending_submission.json',self.pending)
        mb.atomic(self.branch/'ERROR.json',dict(error=str(subprocess.CalledProcessError(1,command)),job='50013914'))
        (self.branch/'logs').mkdir();(self.branch/'logs/control_50013914.err').write_text('sbatch: error: Batch job submission failed: Unexpected message received\n')
        self.account=''.join('50014103_%d|slide_online_f_abc123|test|COMPLETED|0:0|%s\n'%(k,self.branch) for k in range(1,131))
        for k,g in enumerate(self.genes,1):
            (self.td/('worker_50014103_%d.err'%k)).write_text('')
            out=self.d/(g+'.rds');out.write_bytes(b'validated')
            mb.atomic(self.d/'receipts'/(g+'.json'),dict(gene=g,input_sha256=self.task['input_sha256'],files={g+'.rds':dict(stamp=mb.stamp(out))},completed=200))
    def tearDown(self):
        assert (ROOT/'tmp').resolve() in self.root.parents
        self.temp.cleanup()
    @contextlib.contextmanager
    def accounting(self,text=None):
        with patch.dict(os.environ,USER='test'),patch.object(mb,'queue_output',return_value=self.account if text is None else text):yield
    def test_repair_matches_all_130_receipts_and_preserves_files(self):
        before={str(p):mb.sha(p) for p in self.branch.rglob('*') if p.is_file()}
        with self.accounting():r=repair_submission.prepare(self.branch,mb)
        self.assertEqual(r['submission']['job'],'50014103');self.assertEqual(r['workers'],130)
        self.assertEqual(before,{str(p):mb.sha(p) for p in self.branch.rglob('*') if p.is_file()})
    def test_absence_is_not_rejection(self):
        for p in self.td.glob('worker_*.err'):p.unlink()
        with self.accounting():self.assertIsNone(submission_receipts.recover(self.branch,mb,self.pending))
    def test_delayed_accounting_returns_no_ack(self):
        with self.accounting(''):self.assertIsNone(submission_receipts.recover(self.branch,mb,self.pending))
    def test_wrong_identity_is_rejected(self):
        with self.accounting(self.account.replace('|test|','|another|')):
            with self.assertRaisesRegex(RuntimeError,'identity mismatch'):repair_submission.prepare(self.branch,mb)
    def test_changed_result_not_accepted(self):
        (self.d/(self.genes[0]+'.rds')).write_bytes(b'corrupt')
        with self.accounting():
            with self.assertRaisesRegex(RuntimeError,'output changed'):repair_submission.prepare(self.branch,mb)
    def test_multiple_arrays_are_rejected(self):
        (self.td/'worker_50014104_1.err').write_text('')
        with self.accounting():
            with self.assertRaisesRegex(RuntimeError,'Multiple arrays'):repair_submission.prepare(self.branch,mb)
    def test_partial_failure_is_not_cleared(self):
        with self.accounting(self.account.replace('COMPLETED|0:0','FAILED|1:0',1)):
            with self.assertRaisesRegex(RuntimeError,'not completed successfully'):repair_submission.prepare(self.branch,mb)
    def test_submit_lost_ack_recovers_without_second_sbatch(self):
        (self.branch/'pending_submission.json').unlink()
        exc=subprocess.CalledProcessError(1,['sbatch'],stderr='sbatch: error: Batch job submission failed: Unexpected message received')
        with self.accounting(),patch.object(mb.subprocess,'check_output',side_effect=exc) as rpc:
            job=mb.submit(self.branch,'workers',[],self.task)
        self.assertEqual(job,'50014103');self.assertEqual(rpc.call_count,1)
        self.assertFalse((self.branch/'pending_submission.json').exists())
        self.assertEqual(mb.read(self.td/'submission.json')['job'],'50014103')
    def test_bad_resource_request_not_retried(self):
        (self.branch/'pending_submission.json').unlink()
        exc=subprocess.CalledProcessError(1,['sbatch'],stderr='Invalid partition name specified')
        with patch.object(mb.subprocess,'check_output',side_effect=exc) as rpc:
            with self.assertRaises(subprocess.CalledProcessError):mb.submit(self.branch,'workers',[],self.task)
        self.assertEqual(rpc.call_count,1);self.assertTrue((self.branch/'pending_submission.json').exists())
    def test_runtime_latest_success_and_stable_median_fallback(self):
        genes=self.genes[:4]
        mb.atomic(self.d/'attempts'/(genes[0]+'_1.json'),dict(gene=genes[0],stage=93,started=100))
        mb.atomic(self.d/'attempts'/(genes[0]+'_2.json'),dict(gene=genes[0],stage=93,started=170))
        mb.atomic(self.d/'attempts'/(genes[1]+'_3.json'),dict(gene=genes[1],stage=93,started=50))
        rows=runtime_order.collect(self.branch,mb,genes,[dict(directory=str(self.d),stage=93,iteration=23)])
        self.assertEqual([r['seconds'] for r in rows],[30,150,90,90])
        self.assertEqual([r['observed'] for r in rows],[True,True,False,False])
        ordered=runtime_order.order(genes,{r['gene']:r['seconds'] for r in rows})
        self.assertEqual(ordered,[genes[1],genes[2],genes[3],genes[0]])
        self.assertEqual(set(ordered),set(genes))

if __name__=='__main__':unittest.main()
