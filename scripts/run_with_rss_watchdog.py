"""Run a command while polling its macOS resident-memory use.

The child is terminated if its sampled resident set reaches the configured
limit. This is an RSS watchdog, not a kernel allocation cap; polling can miss a
very brief spike between samples. Use for controlled synthetic profiling.
"""

from __future__ import annotations

import argparse
import ctypes
import os
import signal
import subprocess
import sys
import time


class _WatchdogInterrupted(Exception):
    def __init__(self, signum: int):
        self.signum = signum


class _ProcTaskInfo(ctypes.Structure):
    """Subset-compatible 96-byte Darwin ``proc_taskinfo`` layout."""

    _fields_ = [
        ("virtual_size", ctypes.c_uint64),
        ("resident_size", ctypes.c_uint64),
        ("total_user", ctypes.c_uint64),
        ("total_system", ctypes.c_uint64),
        ("threads_user", ctypes.c_uint64),
        ("threads_system", ctypes.c_uint64),
        ("policy", ctypes.c_int32),
        ("faults", ctypes.c_int32),
        ("pageins", ctypes.c_int32),
        ("cow_faults", ctypes.c_int32),
        ("messages_sent", ctypes.c_int32),
        ("messages_received", ctypes.c_int32),
        ("syscalls_mach", ctypes.c_int32),
        ("syscalls_unix", ctypes.c_int32),
        ("context_switches", ctypes.c_int32),
        ("thread_count", ctypes.c_int32),
        ("running_threads", ctypes.c_int32),
        ("priority", ctypes.c_int32),
    ]


def _process_rss_reader():
    library = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
    proc_pidinfo = library.proc_pidinfo
    proc_pidinfo.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_uint64,
        ctypes.POINTER(_ProcTaskInfo),
        ctypes.c_int,
    ]
    proc_pidinfo.restype = ctypes.c_int

    def read_rss(pid: int) -> int:
        info = _ProcTaskInfo()
        result = proc_pidinfo(pid, 4, 0, ctypes.byref(info), ctypes.sizeof(info))
        if result != ctypes.sizeof(info):
            error = ctypes.get_errno()
            raise OSError(error, f"proc_pidinfo failed for child {pid}")
        return int(info.resident_size)

    return read_rss


def _terminate_child_group(child: subprocess.Popen, *, grace_seconds: float = 2.0) -> None:
    """Stop the watched process and any descendants it launched."""
    if child.poll() is not None:
        return
    try:
        os.killpg(child.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    except PermissionError:
        # A just-exited child may disappear between poll() and killpg().
        if child.poll() is None:
            raise
    try:
        child.wait(timeout=grace_seconds)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(child.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        child.wait()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit-gb", type=float, required=True)
    parser.add_argument("--interval-seconds", type=float, default=0.2)
    parser.add_argument("--report-interval-seconds", type=float, default=30.0)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    arguments = parser.parse_args()
    command = arguments.command[1:] if arguments.command[:1] == ["--"] else arguments.command
    if (
        arguments.limit_gb <= 0
        or arguments.interval_seconds <= 0
        or arguments.report_interval_seconds <= 0
        or not command
    ):
        parser.error("provide positive limits, intervals, and a command after --")
    if sys.platform != "darwin":
        parser.error("RSS monitoring currently requires macOS libproc")
    read_rss = _process_rss_reader()
    limit_bytes = int(arguments.limit_gb * 1024**3)
    child = subprocess.Popen(command, start_new_session=True)
    peak_rss = 0
    exceeded = False
    next_report = time.monotonic() + arguments.report_interval_seconds
    previous_handlers = {
        signum: signal.getsignal(signum)
        for signum in (signal.SIGINT, signal.SIGTERM)
    }

    def interrupt(signum: int, _frame) -> None:
        raise _WatchdogInterrupted(signum)

    for signum in previous_handlers:
        signal.signal(signum, interrupt)
    interrupted = None
    monitor_error = None
    try:
        while True:
            if child.poll() is not None:
                break
            try:
                current_rss = read_rss(child.pid)
            except OSError:
                # proc_pidinfo can report ESRCH if the child exits after the
                # loop poll above. In that race the child's exit status is the
                # useful result; only a live child turns this into a monitor
                # failure.
                if child.poll() is not None:
                    break
                raise
            peak_rss = max(peak_rss, current_rss)
            now = time.monotonic()
            if now >= next_report:
                print(
                    f"rss_sample_gb={current_rss / 1024**3:.3f} "
                    f"peak_sampled_rss_gb={peak_rss / 1024**3:.3f} "
                    f"rss_limit_gb={arguments.limit_gb:g}",
                    flush=True,
                )
                next_report = now + arguments.report_interval_seconds
            if current_rss >= limit_bytes:
                exceeded = True
                _terminate_child_group(child)
                break
            time.sleep(arguments.interval_seconds)
    except _WatchdogInterrupted as error:
        interrupted = error.signum
        _terminate_child_group(child)
    except BaseException as error:
        monitor_error = error
        if child.poll() is None:
            _terminate_child_group(child)
    finally:
        for signum, handler in previous_handlers.items():
            signal.signal(signum, handler)
    exit_code = child.wait()
    if monitor_error is not None:
        print(f"rss_monitor_error={monitor_error!r} child_exit={exit_code}", flush=True)
        return 125
    if exceeded:
        exit_code = 124
    print(
        f"child_exit={exit_code} rss_limit_gb={arguments.limit_gb:g} "
        f"peak_sampled_rss_gb={peak_rss / 1024**3:.3f} "
        f"terminated_for_rss_limit={str(exceeded).lower()}",
        flush=True,
    )
    if interrupted is not None:
        print(f"watchdog_interrupted_by_signal={interrupted}", flush=True)
        return 128 + interrupted
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
