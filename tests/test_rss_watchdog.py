from __future__ import annotations

import errno
import importlib.util
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest


@pytest.mark.skipif(sys.platform != "darwin", reason="RSS watchdog uses macOS libproc")
def test_interrupting_watchdog_stops_the_watched_process_group(tmp_path: Path) -> None:
    watchdog = Path(__file__).parents[1] / "scripts" / "run_with_rss_watchdog.py"
    pid_path = tmp_path / "child.pid"
    child_code = (
        "import os, pathlib, time; "
        "pathlib.Path(__import__('sys').argv[1]).write_text(str(os.getpid())); "
        "time.sleep(60)"
    )
    process = subprocess.Popen(
        [
            sys.executable,
            str(watchdog),
            "--limit-gb",
            "16",
            "--interval-seconds",
            "0.02",
            "--",
            sys.executable,
            "-c",
            child_code,
            str(pid_path),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    try:
        deadline = time.monotonic() + 10
        while not pid_path.exists() and time.monotonic() < deadline:
            if process.poll() is not None:
                stdout, stderr = process.communicate()
                pytest.fail(f"watchdog exited early: {stdout}\n{stderr}")
            time.sleep(0.02)
        assert pid_path.exists(), "watched child did not start"
        child_pid = int(pid_path.read_text())

        process.send_signal(signal.SIGINT)
        stdout, stderr = process.communicate(timeout=10)

        assert process.returncode == 128 + signal.SIGINT, stdout + stderr
        assert "watchdog_interrupted_by_signal=" in stdout
        with pytest.raises(ProcessLookupError):
            os.kill(child_pid, 0)
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)


@pytest.mark.skipif(sys.platform != "darwin", reason="RSS watchdog uses macOS libproc")
def test_child_exit_between_poll_and_rss_read_is_reported_as_normal_exit(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    watchdog = Path(__file__).parents[1] / "scripts" / "run_with_rss_watchdog.py"
    spec = importlib.util.spec_from_file_location("rss_watchdog_race_test", watchdog)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    real_popen = subprocess.Popen

    class PollRaceProxy:
        def __init__(self, *args, **kwargs) -> None:
            self.child = real_popen(*args, **kwargs)
            self.first_poll = True

        def poll(self):
            if self.first_poll:
                self.first_poll = False
                return None
            return self.child.poll()

        def __getattr__(self, name):
            return getattr(self.child, name)

    def exited_during_read(_pid: int) -> int:
        time.sleep(0.05)
        raise ProcessLookupError(errno.ESRCH, "child exited before RSS read")

    monkeypatch.setattr(module.subprocess, "Popen", PollRaceProxy)
    monkeypatch.setattr(module, "_process_rss_reader", lambda: exited_during_read)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(watchdog), "--limit-gb", "16", "--interval-seconds", "0.01",
            "--report-interval-seconds", "1", "--", sys.executable, "-c",
            "__import__('time').sleep(0.3)",
        ],
    )

    assert module.main() == 0
    output = capsys.readouterr().out
    assert "child_exit=0" in output
    assert "rss_monitor_error=" not in output


@pytest.mark.skipif(sys.platform != "darwin", reason="RSS watchdog uses macOS libproc")
def test_transient_libproc_lookup_failure_does_not_kill_a_live_child(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    watchdog = Path(__file__).parents[1] / "scripts" / "run_with_rss_watchdog.py"
    spec = importlib.util.spec_from_file_location("rss_watchdog_transient_test", watchdog)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    reads = 0

    def transiently_missing(_pid: int) -> int:
        nonlocal reads
        reads += 1
        if reads == 1:
            raise ProcessLookupError(errno.ESRCH, "transient proc_pidinfo lookup failure")
        return 1024 * 1024

    monkeypatch.setattr(module, "_process_rss_reader", lambda: transiently_missing)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(watchdog), "--limit-gb", "16", "--interval-seconds", "0.01",
            "--report-interval-seconds", "10", "--", sys.executable, "-c",
            "__import__('time').sleep(1.3)",
        ],
    )

    assert module.main() == 0
    output = capsys.readouterr().out
    assert reads > 1
    assert "child_exit=0" in output
    assert "rss_monitor_error=" not in output
