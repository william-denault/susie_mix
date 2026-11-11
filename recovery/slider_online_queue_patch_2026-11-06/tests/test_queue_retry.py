import contextlib, errno, importlib.util, io, os, subprocess, sys, unittest
from pathlib import Path
from unittest.mock import patch

HERE=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('queue_runtime',HERE/'runtime/minibatch.py')
mb=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mb;spec.loader.exec_module(mb)

class QueueTests(unittest.TestCase):
    def setUp(self):
        self.command=['squeue','--noheader','--array','--states=all','--user','test',
                      '--format=%i|%F|%100j|%T|%r|%S']
        self.row='123_1|123|slide_online_f_abc123|RUNNING|node|N/A\n'
        self.timeout=subprocess.TimeoutExpired(self.command,60,output='partial, untrusted output')

    def test_timeout_retries_actual_queue_read_not_empty_snapshot(self):
        with patch.dict(os.environ,USER='test'),patch.object(mb.subprocess,'check_output',
             side_effect=[self.timeout,self.timeout,self.row]) as run,patch.object(mb.time,'sleep') as sleep, \
             contextlib.redirect_stderr(io.StringIO()) as log:
            rows=mb.queue()
        self.assertEqual([r['job'] for r in rows],['123_1'])
        self.assertEqual([c.args[0] for c in sleep.call_args_list],[10,20])
        self.assertEqual(run.call_count,3)
        self.assertTrue(all(c.kwargs['timeout']==60 for c in run.call_args_list))
        self.assertIn('recovered after 3 attempts',log.getvalue())

    def test_known_transient_failures_retry(self):
        cases=[subprocess.CalledProcessError(1,self.command,stderr='slurm_load_jobs error: Socket timed out on send/recv operation'),
               subprocess.CalledProcessError(1,self.command,output=b'Unable to contact slurm controller (connect failure)'),
               OSError(errno.EAGAIN,'Resource temporarily unavailable')]
        for error in cases:
            with self.subTest(error=error),patch.dict(os.environ,USER='test'), \
                 patch.object(mb.subprocess,'check_output',side_effect=[error,self.row]), \
                 patch.object(mb.time,'sleep'),contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(len(mb.queue()),1)

    def test_permanent_or_malformed_response_stops(self):
        cases=[subprocess.CalledProcessError(1,self.command,stderr='squeue: error: Invalid user: test'),
               FileNotFoundError(errno.ENOENT,'squeue missing'), 'garbled response']
        for error in cases:
            with self.subTest(error=error),patch.dict(os.environ,USER='test'), \
                 patch.object(mb.subprocess,'check_output',side_effect=[error] if isinstance(error,Exception) else None,
                              return_value=error),patch.object(mb.time,'sleep') as sleep:
                with self.assertRaises(Exception):mb.queue()
                sleep.assert_not_called()

    def test_controller_waits_through_long_outage_at_same_operation(self):
        outcomes=[self.timeout]*12+[self.row]
        with patch.dict(os.environ,USER='test'),patch.object(mb.subprocess,'check_output',side_effect=outcomes), \
             patch.object(mb.time,'sleep') as sleep,patch.object(mb.time,'monotonic',side_effect=range(0,200000,1000)), \
             contextlib.redirect_stderr(io.StringIO()),mb.wait_for_scheduler():
            self.assertEqual(len(mb.queue()),1)
        self.assertEqual(sleep.call_count,12)
        self.assertEqual([c.args[0] for c in sleep.call_args_list], [10,20,40]+[60]*9)
        self.assertFalse(mb.QUEUE_WAIT_FOREVER)

    def test_interactive_outage_has_bounded_wait(self):
        with patch.dict(os.environ,USER='test'),patch.object(mb.subprocess,'check_output',side_effect=self.timeout), \
             patch.object(mb.time,'sleep'),patch.object(mb.time,'monotonic',side_effect=[0,0,301]), \
             contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaisesRegex(RuntimeError,'no queue snapshot'):mb.queue()

    def test_control_enables_wait_and_restores_context_even_on_error(self):
        def action(branch):
            self.assertTrue(mb.QUEUE_WAIT_FOREVER)
            raise RuntimeError('actual validation error')
        with patch.object(mb,'control_once',side_effect=action):
            with self.assertRaisesRegex(RuntimeError,'validation error'):mb.control(Path('unused'))
        self.assertFalse(mb.QUEUE_WAIT_FOREVER)

if __name__=='__main__':unittest.main()
