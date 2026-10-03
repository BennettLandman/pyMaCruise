from __future__ import annotations

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
