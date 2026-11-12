# bactipipe/scripts/runner.py
from __future__ import annotations
from pathlib import Path
from typing import List, Optional, Tuple
import logging
import math
import os
import shutil
import signal
import subprocess
import threading


def _bounded_run(cmd, cwd, timeout, grace=5):
    """Bound a tool and its descendants, including conda/KMA grandchildren.

    Invoked by the sequential typing driver on its main thread. Forward parent
    cancellation too: a separate process group must not outlive a stopped cohort.
    """
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("Tool timeout must be a finite positive number")
    if threading.current_thread() is not threading.main_thread():
        raise ValueError("Bounded tools must run on the main thread for cancellation safety")
    previous = {}
    proc = None

    def cancel(signum, frame):
        raise SystemExit(128 + signum)

    def kill_group(sig):
        try:
            os.killpg(proc.pid, sig)
        except ProcessLookupError:
            pass

    def terminate():
        # Ignore repeat cancellation while reaping the entire tool process group.
        for sig in previous:
            signal.signal(sig, signal.SIG_IGN)
        kill_group(signal.SIGTERM)
        try:
            proc.communicate(timeout=grace)
        except subprocess.TimeoutExpired:
            pass
        finally:
            # The group leader can exit before a grandchild which ignored TERM.
            kill_group(signal.SIGKILL)
        return proc.communicate()[0]

    try:
        for sig in (signal.SIGTERM, signal.SIGINT):
            previous[sig] = signal.signal(sig, cancel)
        proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True,
                                start_new_session=True)
        try:
            output, _ = proc.communicate(timeout=timeout)
            return proc.returncode, output
        except subprocess.TimeoutExpired:
            output = terminate()
            return 124, (output or "") + f"\n[runner] Timed out after {timeout:g} seconds; tool process group terminated.\n"
        except BaseException:
            terminate()
            raise
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)

def compose_tool_cmd(cmd: List[str], env_name: Optional[str]) -> List[str]:
    """Wrap command to run inside a conda/mamba/micromamba env if requested."""
    if env_name:
        if shutil.which("conda"):
            return ["conda", "run", "--no-capture-output", "-n", env_name] + cmd
        if shutil.which("micromamba"):
            return ["micromamba", "run", "-n", env_name] + cmd
        if shutil.which("mamba"):
            return ["mamba", "run", "-n", env_name] + cmd
    return cmd

def run(cmd, cwd, log, env_name=None, timeout=None):
    cmd_exec = compose_tool_cmd(cmd, env_name)
    # normalize everything to str for safe logging/subprocess
    if any(isinstance(x, bool) for x in cmd_exec):
        log.error(f"BUG: command contains boolean(s): {cmd_exec!r}")
    cmd_exec = [str(x) for x in cmd_exec]
    log.debug(f"RUN: {' '.join(cmd_exec)}")
    try:
        if timeout is not None:
            return _bounded_run(cmd_exec, str(cwd) if cwd else None, timeout)
        p = subprocess.run(
            cmd_exec,
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        return p.returncode, p.stdout
    except Exception as e:
        return 127, f"[runner] Failed to execute: {e}"

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)
