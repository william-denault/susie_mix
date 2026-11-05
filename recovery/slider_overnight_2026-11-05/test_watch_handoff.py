import os
from pathlib import Path
import subprocess
import tempfile
import unittest

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
SH = Path('C:/Users/willi/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/git/usr/bin/sh.exe')


class WatchTests(unittest.TestCase):
    def run_case(self, messages, expected_calls, expected_exit, marker=None):
        scratch = ROOT / 'tmp'
        scratch.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='overnight_watch_', dir=str(scratch)) as tmp:
            directory = Path(tmp).resolve()
            self.assertIn(scratch.resolve(), directory.parents)
            bundle = directory / 'bundle'
            branch = directory / 'branch'
            bundle.mkdir()
            branch.mkdir()
            if marker:
                (branch / marker).write_bytes(b'preserve-me')
            lines = ['#!/bin/bash', 'n=0', '[ ! -f calls.txt ] || read -r n < calls.txt',
                     'n=$((n + 1))', 'echo "$n" > calls.txt', 'case "$n" in']
            for number, (message, code) in enumerate(messages, 1):
                quoted = "'" + message.replace("'", "'\"'\"'") + "'"
                lines.append('{} ) printf "%s\\n" {}; exit {} ;;'.format(number, quoted, code))
            lines += ['*) echo UNEXPECTED_EXTRA_LAUNCH; exit 77 ;;', 'esac']
            (bundle / 'start.sh').write_bytes(('\n'.join(lines) + '\n').encode())
            env = os.environ.copy()
            env['PATH'] = str(SH.parent) + os.pathsep + env.get('PATH', '')
            env['SLIDER_WATCH_RETRY_SECONDS'] = '0'
            env['WATCH_TEST_SH'] = SH.as_posix()
            command = 'bash() { "$WATCH_TEST_SH" "$@"; }; export -f bash; source "$1" "$2" "$3"'
            result = subprocess.run([str(SH), '-c', command, 'watch-test',
                                     (BASE / 'watch_handoff.sh').as_posix(),
                                     bundle.as_posix(), branch.as_posix()],
                                    env=env, text=True, capture_output=True, timeout=20)
            self.assertEqual(result.returncode, expected_exit, result.stdout + result.stderr)
            calls = bundle / 'calls.txt'
            self.assertEqual(int(calls.read_text()) if calls.exists() else 0, expected_calls,
                             result.stdout + result.stderr)
            if marker:
                self.assertEqual((branch / marker).read_bytes(), b'preserve-me')
            return result.stdout

    def test_waits_and_existing_launcher_lock_then_success(self):
        messages = [(message, 1) for message in (
            'STOP: [Errno 11] Resource temporarily unavailable',
            'STOP: No safe pending-controller window; retry --resume --handoff',
            'STOP: No safe submission window/free slots; retry --resume later')]
        messages.append(('Submitted controller 12345', 0))
        self.assertIn('SUCCESS:', self.run_case(messages, 4, 0))

    def test_existing_controller_exits_without_a_second_attempt(self):
        self.run_case([('STOP: A branch controller is already queued/running', 1)], 1, 0)

    def test_stops_on_real_errors(self):
        for error in ('Cannot take over 014: later iteration folders already exist: iteration_014',
                      'Bundle file changed: engine.R',
                      'Inspect ambiguous pending_submission.json before retrying',
                      'Recorded controller is no longer safely pending; no cancellation performed'):
            with self.subTest(error=error):
                self.run_case([('STOP: ' + error, 1)], 1, 1)

    def test_existing_error_is_preserved_without_resume(self):
        self.run_case([], 0, 1, 'ERROR.json')

    def test_ambiguous_submission_not_masked_by_existing_job(self):
        self.run_case([('STOP: A branch controller is already queued/running', 1)],
                      1, 1, 'pending_submission.json')


if __name__ == '__main__':
    unittest.main()
