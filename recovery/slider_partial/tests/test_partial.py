import hashlib,json,os,sys,tempfile,types,unittest
from pathlib import Path
from unittest import mock
sys.dont_write_bytecode=True;sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import storage,partial_resume as partial,handoff

class PartialTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir='tmp');self.root=Path(self.tmp.name).resolve()
        self.project=self.root/'project';self.project.mkdir();self.scratch=self.root/'scratch';self.scratch.mkdir()
        self.rel='results_slide_prior_em/iteration_004/results/GENE.rds';self.src=self.project/self.rel
        self.src.parent.mkdir(parents=True);self.old=b'fit data\n'+b'redundant call '*300;self.src.write_bytes(self.old)
        self.stage=self.scratch/hashlib.sha256(self.rel.encode()).hexdigest()[:24]
        self.ledger=storage.Ledger(self.root/'receipts.sqlite')
    def tearDown(self):self.ledger.db.close();self.tmp.cleanup()
    def compact(self,src,dst,report):dst.write_bytes(b'fit data\n');report.write_text('removed_calls\n1\n')
    def recover(self):return partial.recover_transactions(self.project,self.scratch,self.ledger,{self.rel})
    def test_start_retains_original_without_new_compaction(self):
        self.stage.mkdir();t=dict(relative=self.rel,stage='start',old_hash=storage.hashes(self.src),old_bytes=len(self.old),original_stamp=storage.stamp(self.src))
        storage.atomic_json(self.stage/'transaction.json',t);(self.stage/'compact.rds').write_bytes(b'incomplete')
        self.recover();self.assertEqual(self.src.read_bytes(),self.old)
        self.assertEqual(self.ledger.get(self.rel)['removed_calls'],0)
    def test_vacated_transaction_finishes_existing_verified_compact(self):
        real=storage.copy_verified
        def fail(src,dst,expected=None):
            if Path(dst).name.startswith('.storage_compact_'):raise OSError('interrupted')
            return real(src,dst,expected)
        with mock.patch.object(storage,'copy_verified',side_effect=fail):
            with self.assertRaises(OSError):storage.transact(self.project,self.rel,self.stage,self.ledger,self.compact)
        self.assertFalse(self.src.exists());self.recover();self.assertEqual(self.src.read_bytes(),b'fit data\n')
    def test_replaced_but_uncommitted_is_adopted(self):
        with mock.patch.object(self.ledger,'put',side_effect=OSError('interrupted')):
            with self.assertRaises(OSError):storage.transact(self.project,self.rel,self.stage,self.ledger,self.compact)
        self.recover();self.assertEqual(self.src.read_bytes(),b'fit data\n');self.assertIsNotNone(self.ledger.get(self.rel))
    def test_untouched_files_are_not_read_or_recompressed(self):
        with mock.patch.object(storage,'hashes',side_effect=AssertionError('must not read untouched fit')):
            self.assertEqual(self.recover(),0);self.assertEqual(partial.verify_committed(self.project,self.ledger,{self.rel}),[])
        self.assertEqual(self.src.read_bytes(),self.old)
    def test_changed_start_original_blocks_recovery(self):
        self.stage.mkdir();t=dict(relative=self.rel,stage='start',old_hash=storage.hashes(self.src),old_bytes=len(self.old),original_stamp=storage.stamp(self.src))
        storage.atomic_json(self.stage/'transaction.json',t);self.src.write_bytes(b'external edit')
        with self.assertRaisesRegex(RuntimeError,'original changed'):self.recover()
        self.assertTrue((self.stage/'transaction.json').exists())
    def test_changed_committed_file_is_rejected(self):
        storage.transact(self.project,self.rel,self.stage,self.ledger,self.compact);self.src.write_bytes(b'changed')
        with self.assertRaisesRegex(RuntimeError,'Committed result changed'):partial.verify_committed(self.project,self.ledger,{self.rel})
    def test_corrupt_escrow_blocks_finish(self):
        real=storage.copy_verified
        def fail(src,dst,expected=None):
            if Path(dst).name.startswith('.storage_compact_'):raise OSError('interrupted')
            return real(src,dst,expected)
        with mock.patch.object(storage,'copy_verified',side_effect=fail):
            with self.assertRaises(OSError):storage.transact(self.project,self.rel,self.stage,self.ledger,self.compact)
        (self.stage/'compact.rds').write_bytes(b'corrupt')
        with self.assertRaisesRegex(RuntimeError,'Escrow'):self.recover()
        self.assertTrue((self.stage/'original.rds').exists())
    def test_inflight_path_must_belong_to_original_plan(self):
        self.stage.mkdir();storage.atomic_json(self.stage/'transaction.json',dict(relative='outside',stage='start'))
        with self.assertRaisesRegex(RuntimeError,'Unexpected scratch'):self.recover()

class HandoffTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir='tmp');self.root=Path(self.tmp.name).resolve()
        self.original=self.root/'original';self.original.mkdir();(self.original/'queued_jobs.tsv').write_text('compact\t85\nresume\t86\n')
        self.rows=[['85','slide_storage_compact','RUNNING'],['86','slide_storage_resume','PENDING'],['11','simulation','RUNNING']]
    def tearDown(self):self.tmp.cleanup()
    def test_only_recorded_jobs_selected_in_pending_resume_first_order(self):
        _,stop=handoff.identify(self.original,self.rows);self.assertEqual(stop,['86','85'])
    def test_running_resume_or_other_em_refuses_cancellation(self):
        for rows in (self.rows+[['12','slide_parallel_fit','RUNNING']], [self.rows[0],['86','slide_storage_resume','RUNNING']]):
            with self.assertRaises(RuntimeError):handoff.identify(self.original,rows)
    def test_wrong_job_name_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError,'identity differs'):handoff.identify(self.original,[['85','unrelated','RUNNING']])
    def test_held_job_is_created_before_scoped_cancellations_then_released(self):
        project=self.root/'project';project.mkdir();scratch=self.root/'scratch';scratch.mkdir();state=self.root/'state'
        storage.atomic_json(self.original/'identity.json',dict(project=str(project),kind='compact',expected=5,scratch=str(scratch)))
        storage.atomic_json(scratch/'owner.json',dict(project=str(project),state=str(self.original)))
        for name in ('snapshot.rds','files.json','receipts.sqlite'):(self.original/name).write_text('fixture')
        args=types.SimpleNamespace(project=project,original_state=self.original,scratch=scratch,state=state)
        fake=types.SimpleNamespace(LOCK_EX=1,LOCK_NB=2,flock=lambda *a:None)
        events=[]
        def submit(cmd,**kw):events.append(cmd);return '99\n'
        def call(cmd,**kw):events.append(cmd);return 0
        with mock.patch.dict(sys.modules,{'fcntl':fake}),mock.patch.object(storage,'verify_bundle'),mock.patch.object(storage,'check_sources'),\
             mock.patch.object(storage,'queue',side_effect=[self.rows,self.rows+[['99','slide_partial_finalize','PENDING']]]),\
             mock.patch.object(handoff.subprocess,'check_output',side_effect=submit),mock.patch.object(handoff.subprocess,'check_call',side_effect=call):
            handoff.run(args)
        self.assertIn('--hold',events[0]);self.assertIn('--dependency=afterany:86:85',events[0])
        self.assertEqual(events[1:], [['scancel','--state=PENDING','86'],['scancel','85'],['scontrol','release','99']])
    def test_live_maintenance_blocks_finalization(self):
        self.assertTrue(partial.unsafe_jobs([['85','slide_storage_compact','RUNNING']],'99'))
        self.assertFalse(partial.unsafe_jobs([['99','slide_partial_finalize','RUNNING'],['11','simulation','RUNNING']],'99'))

if __name__=='__main__':unittest.main()
