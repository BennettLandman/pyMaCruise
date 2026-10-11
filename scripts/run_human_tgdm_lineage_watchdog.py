#!/usr/bin/env python3
"""Run the saved-crop TGDM lineage diagnostic with RSS and time bounds."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-attempt-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--limit-gb", type=float, default=8.0)
    parser.add_argument("--timeout-seconds", type=float, default=600.0)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.limit_gb <= 0 or args.timeout_seconds <= 0:
        parser.error("RSS and timeout limits must be positive")

    repo_root = Path(__file__).resolve().parents[1]
    watchdog = Path(__file__).with_name("run_with_rss_watchdog.py")
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()) and not args.overwrite:
        parser.error(f"output directory is not empty: {output_dir}; use --overwrite")
    log_path = output_dir.with_name(output_dir.name + ".watchdog.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        str(watchdog),
        "--limit-gb",
        str(args.limit_gb),
        "--interval-seconds",
        "0.2",
        "--report-interval-seconds",
        "10",
        "--",
        sys.executable,
        "-m",
        "pymacruise.tgdm_diagnostic",
        "--core-attempt-dir",
        str(args.core_attempt_dir.resolve()),
        "--output-dir",
        str(output_dir),
        *( ["--overwrite"] if args.overwrite else [] ),
    ]
    environment = os.environ.copy()
    environment.setdefault("NUMBA_NUM_THREADS", "1")
    started = time.monotonic()
    timed_out = False
    with log_path.open("w", encoding="utf-8") as log:
        log.write("command=" + json.dumps(command) + "\n")
        log.flush()
        process = subprocess.Popen(
            command,
            cwd=repo_root,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            return_code = process.wait(timeout=args.timeout_seconds)
        except subprocess.TimeoutExpired:
            timed_out = True
            process.send_signal(signal.SIGINT)
            try:
                return_code = process.wait(timeout=20.0)
            except subprocess.TimeoutExpired:
                process.kill()
                return_code = process.wait()
    runtime = time.monotonic() - started
    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    sampled = re.findall(r"peak_sampled_rss_gb=([0-9.]+)", log_text)
    rss_samples = re.findall(r"rss_sample_gb=([0-9.]+)", log_text)
    watchdog_match = re.search(
        r"terminated_for_rss_limit=(true|false)", log_text
    )
    watchdog_summary = {
        "limit_gib": args.limit_gb,
        "timeout_seconds": args.timeout_seconds,
        "runtime_seconds": runtime,
        "peak_sampled_rss_gib": float(sampled[-1]) if sampled else None,
        "last_sampled_rss_gib": float(rss_samples[-1]) if rss_samples else None,
        "terminated_for_rss_limit": (
            watchdog_match.group(1) == "true" if watchdog_match else None
        ),
        "timed_out": timed_out,
        "child_exit_code": return_code,
        "sampling_interval_seconds": 0.2,
        "watchdog_log": str(log_path),
        "limitation": "sampled RSS watchdog can miss a brief peak between samples",
    }

    report_path = output_dir / "lineage.json"
    if report_path.is_file():
        report = json.loads(report_path.read_text(encoding="utf-8"))
        report["watchdog"] = watchdog_summary
        report["runtime_seconds_including_watchdog"] = runtime
        report_path.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(report_path)
    else:
        output_dir.mkdir(parents=True, exist_ok=True)
        watchdog_summary["status"] = (
            "timed_out" if timed_out else
            "rss_limit_exceeded" if return_code == 124 else
            "diagnostic_failed_or_interrupted"
        )
        (output_dir / "watchdog.json").write_text(
            json.dumps(watchdog_summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(output_dir / "watchdog.json")
    print(
        " ".join(
            (
                f"exit={return_code}",
                f"rss_peak_gib={watchdog_summary['peak_sampled_rss_gib']}",
                f"runtime_seconds={runtime:.1f}",
                f"timed_out={str(timed_out).lower()}",
            )
        ),
        flush=True,
    )
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
