import logging
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from bactipipe.scripts.runner import _bounded_run, run


def test_bounded_success_and_error():
    assert _bounded_run([sys.executable, '-c', 'print("ok")'], None, 5) == (0, 'ok\n')
    rc, output = _bounded_run([sys.executable, '-c', 'print("bad");exit(3)'], None, 5)
    assert (rc, output) == (3, 'bad\n')


def test_timeout_preserves_output_and_restores_handlers():
    before = {sig: signal.getsignal(sig) for sig in (signal.SIGTERM, signal.SIGINT)}
    rc, output = run([sys.executable, '-u', '-c', 'import time;print("started");time.sleep(60)'],
                     None, logging.getLogger(__name__), timeout=0.15)
    assert rc == 124 and 'started' in output and 'Timed out' in output
    assert all(signal.getsignal(sig) == value for sig, value in before.items())


def _alive(pid):
    path = Path(f'/proc/{pid}/stat')
    return path.exists() and path.read_text().split(') ', 1)[1][0] != 'Z'


def test_timeout_kills_term_ignoring_grandchild(tmp_path):
    pidfile = tmp_path / 'child.pid'
    child = f'import os,signal,time;signal.signal(signal.SIGTERM,signal.SIG_IGN);open({str(pidfile)!r},"w").write(str(os.getpid()));time.sleep(60)'
    parent = f'import subprocess,time;subprocess.Popen([{sys.executable!r},"-c",{child!r}]);time.sleep(60)'
    rc, _ = _bounded_run([sys.executable, '-c', parent], None, .5, grace=.1)
    assert rc == 124
    assert not _alive(int(pidfile.read_text()))


def test_parent_cancellation_kills_tool(tmp_path):
    pidfile = tmp_path / 'tool.pid'
    child = f'import os,time;open({str(pidfile)!r},"w").write(str(os.getpid()));time.sleep(60)'
    driver = f'from bactipipe.scripts.runner import _bounded_run;_bounded_run([{sys.executable!r},"-c",{child!r}],None,30)'
    parent = subprocess.Popen([sys.executable, '-c', driver])
    try:
        deadline = time.monotonic() + 5
        while not pidfile.exists() and time.monotonic() < deadline:
            time.sleep(.02)
        assert pidfile.exists()
        parent.send_signal(signal.SIGTERM)
        assert parent.wait(timeout=10) == 143
        assert not _alive(int(pidfile.read_text()))
    finally:
        if parent.poll() is None:
            parent.kill()
            parent.wait()


def test_invalid_timeout_does_not_start_tool():
    rc, message = run([sys.executable, '-c', 'raise AssertionError()'], None,
                      logging.getLogger(__name__), timeout=0)
    assert rc == 127 and 'positive' in message
