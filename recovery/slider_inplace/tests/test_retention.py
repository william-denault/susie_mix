import csv,importlib.util,json,os,sys,tempfile,unittest
from pathlib import Path
from unittest import mock
from types import SimpleNamespace
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import retention,storage,launch,orchestrate,common

class RetirementTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir=Path.cwd()/'tmp');self.root=Path(self.tmp.name).resolve()
        assert self.root.parent==(Path.cwd()/'tmp').resolve()
        self.branch=self.root/'branch';self.branch.mkdir();self.source=self.root/'source';self.source.mkdir()
        storage.atomic_json(self.branch/'identity.json',dict(source=str(self.source)))
        self.state=self.branch/'state';self.state.mkdir()
        self.d=self.branch/'results_slide_prior_em/iteration_005';(self.d/'results').mkdir(parents=True)
        self.keep=self.branch/'results_slide_prior_em/iteration_006/results';self.keep.mkdir(parents=True)
        (self.keep/'KEEP.rds').write_bytes(b'keep')
        (self.d/'results/A.rds').write_bytes(b'old A');(self.d/'results/B.rds').write_bytes(b'old B')
        (self.state/'snapshot.rds').write_bytes(b'metadata')
        dirs=[]
        for i in range(6):
            d=self.branch/'results_slide_prior_em'/('iteration_%03d'%i)
            (d/'results').mkdir(parents=True,exist_ok=True);dirs.append(d)
        with (self.state/'directories.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=['directory']);w.writeheader()
            for d in dirs:w.writerow(dict(directory=str(d)))
    def tearDown(self): self.tmp.cleanup()
    def plan(self): return retention.plan_files(self.branch,self.state,7)
    def run_plan(self,**kw): return retention.execute_plan(self.branch,self.state,7,**kw)
    def test_interrupted_deletion_resumes_and_keeps_previous(self):
        self.plan()
        def fail(p): raise RuntimeError('interrupt')
        with self.assertRaisesRegex(RuntimeError,'interrupt'): self.run_plan(after_unlink=fail)
        self.run_plan();self.run_plan()
        self.assertEqual(list((self.d/'results').iterdir()),[])
        self.assertEqual((self.keep/'KEEP.rds').read_bytes(),b'keep')
    def test_replacement_or_modified_old_file_blocks_all_deletion(self):
        self.plan();(self.d/'results/B.rds').write_bytes(b'changed')
        with self.assertRaisesRegex(RuntimeError,'changed'): self.run_plan()
        self.assertTrue((self.d/'results/A.rds').exists())
    def test_unexpected_file_blocks_all_deletion(self):
        self.plan();(self.d/'results/unexpected.txt').write_text('keep')
        with self.assertRaisesRegex(RuntimeError,'Unplanned'): self.run_plan()
        self.assertTrue((self.d/'results/A.rds').exists())
    def test_metadata_change_blocks_deletion(self):
        self.plan();(self.state/'snapshot.rds').write_bytes(b'changed')
        with self.assertRaisesRegex(RuntimeError,'metadata changed'): self.run_plan()
        self.assertTrue((self.d/'results/A.rds').exists())
    def test_current_directory_is_never_accepted(self):
        p=self.plan();p['directories']=[str(self.keep.parent)];storage.atomic_json(self.state/'plan.json',p)
        with self.assertRaisesRegex(RuntimeError,'another run'): self.run_plan()
        self.assertTrue((self.keep/'KEEP.rds').exists())
    def test_redirected_result_directory_refused(self):
        self.plan()
        with mock.patch.object(Path,'is_symlink',lambda p:str(p)==str(self.d/'results')):
            with self.assertRaisesRegex(RuntimeError,'redirected'): self.run_plan()
    def test_prior_and_pooled_files_are_not_retired(self):
        (self.d/'priors.csv').write_bytes(b'prior');(self.d/'recovery/parallel').mkdir(parents=True)
        (self.d/'recovery/parallel/pooled.rds').write_bytes(b'counts')
        (self.d/'recovery/parallel/COMPLETE.rds').write_bytes(b'cert')
        p=self.plan();self.run_plan();retention.archive_certificates(p);retention.archive_certificates(p)
        self.assertEqual((self.d/'priors.csv').read_bytes(),b'prior')
        self.assertEqual((self.d/'recovery/parallel/pooled.rds').read_bytes(),b'counts')
        self.assertEqual((self.d/'recovery/parallel/COMPLETE.before_rolling_retirement.rds').read_bytes(),b'cert')
    def test_live_workers_or_original_jobs_block_cleanup(self):
        with mock.patch.dict(os.environ,{'SLURM_JOB_ID':'123'}):
            own=['123','slide_inplace_control','RUNNING']
            for name in ('slide_expanded_fit','slide_parallel_control','slide_storage_compact'):
                with mock.patch.object(storage,'queue',return_value=[own,['124',name,'PENDING']]):
                    with self.assertRaisesRegex(RuntimeError,'still active'): retention.check_queue()
            with mock.patch.object(storage,'queue',return_value=[own,['125','simulation','RUNNING']]): retention.check_queue()
    def test_gate_waits_for_current_jobs_and_refuses_duplicates(self):
        rows=[['1','slide_storage_resume','PENDING'],['2','slide_parallel_fit','RUNNING'],['3','simulation','RUNNING']]
        self.assertEqual(len(launch.classify(rows)),2)
        for name in ('slide_next_gate','slide_next_setup','slide_expanded_control','slide_inplace_control'):
            with self.assertRaisesRegex(RuntimeError,'duplicate'): launch.classify([['4',name,'PENDING']])
    def test_final_iteration_keeps_039_and_040(self):
        dirs,_=retention.allowed_dirs(self.branch,40)
        self.assertEqual(dirs[-1].name,'iteration_038')
        self.assertFalse(any(d.name in ('iteration_039','iteration_040') for d in dirs))
    def test_partial_finalizer_and_continuation_gate_do_not_block_each_other(self):
        self.assertEqual(common.unsafe([['1','slide_next_gate','PENDING']]),[])
        p=Path(__file__).resolve().parents[1]/'partial_recovery/partial_resume.py'
        spec=importlib.util.spec_from_file_location('partial_test',p);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        self.assertEqual(module.unsafe_jobs([['1','slide_next_gate','PENDING']]),[])
    def test_partial_and_full_recovery_choice(self):
        original=self.root/'original';partial=self.root/'partial';original.mkdir();partial.mkdir()
        self.assertEqual(orchestrate.recovery_choice(original,partial),(partial,True))
        storage.atomic_json(original/'COMPLETE.json',{})
        self.assertEqual(orchestrate.recovery_choice(original,partial),(original,False))
        storage.atomic_json(partial/'handoff.json',dict(phase='held'))
        with self.assertRaisesRegex(RuntimeError,'held'):orchestrate.recovery_choice(original,partial)
        storage.atomic_json(partial/'handoff.json',dict(phase='released'))
        self.assertEqual(orchestrate.recovery_choice(original,partial),(partial,False))
    def test_only_finalized_metadata_can_start_continuation(self):
        state=self.root/'recovery';state.mkdir()
        receipt=dict(project=str(self.branch),latest=5,saved_bytes=123,mode='partial_compaction_finalized',metadata_finalized=False)
        storage.atomic_json(state/'COMPLETE.json',receipt)
        with self.assertRaisesRegex(RuntimeError,'unfinished'):common.recovery_ready(self.branch,state)
        receipt['metadata_finalized']=True;storage.atomic_json(state/'COMPLETE.json',receipt)
        common.recovery_ready(self.branch,state)
    def test_external_dependency_prevents_retirement(self):
        storage.atomic_json(self.branch/'results_slide_prior_em/storage_dependencies.json',[dict(branch='somewhere')])
        with self.assertRaisesRegex(RuntimeError,'depends'):common.check_dependencies(self.branch)
    def test_first_boundary_deletes_only_through_004(self):
        dirs,_=retention.allowed_dirs(self.branch,6)
        self.assertEqual([p.name for p in dirs],['iteration_%03d'%i for i in range(5)])
    def test_gate_reserves_successor_with_shared_lock_and_job_limit(self):
        fake=SimpleNamespace(LOCK_EX=1,LOCK_NB=2,LOCK_UN=4,flock=mock.Mock())
        full=[[str(i),'simulation','RUNNING'] for i in range(300)]
        with mock.patch.dict(sys.modules,{'fcntl':fake}),mock.patch.object(storage,'queue',side_effect=[full,full[:-1]]),mock.patch.object(launch.time,'sleep') as sleep:
            with launch.submission_slot(self.branch):pass
            self.assertEqual(sleep.call_count,1)
            self.assertEqual(fake.flock.call_count,4)
    def test_gate_does_not_submit_above_300(self):
        fake=SimpleNamespace(LOCK_EX=1,LOCK_NB=2,LOCK_UN=4,flock=mock.Mock())
        full=[[str(i),'simulation','RUNNING'] for i in range(300)]
        with mock.patch.dict(sys.modules,{'fcntl':fake}),mock.patch.object(storage,'queue',return_value=full),mock.patch.object(launch.time,'sleep'):
            with self.assertRaisesRegex(RuntimeError,'No free submission slot'):
                with launch.submission_slot(self.branch):self.fail('No slot should be reserved')

if __name__=='__main__':unittest.main()
