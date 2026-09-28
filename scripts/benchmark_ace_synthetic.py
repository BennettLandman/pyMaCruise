"""Measure ACE's default 25-iteration path on deterministic synthetic arrays.

Run from the project root with ``.venv/bin/python scripts/benchmark_ace_synthetic.py``.
Each size is a fresh child process so its peak-RSS measurement is independent.
No input image or benchmark output is written to disk.
"""

from __future__ import annotations

import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def _peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    # macOS reports bytes; Linux and the BSDs used by this project report KiB.
    return value if sys.platform == "darwin" else value * 1024


def _run_size(size: int) -> dict[str, object]:
    from pymacruise.ace import ace_enhance_arrays

    coordinate = np.arange(size, dtype=np.float32) - np.float32((size - 1) / 2)
    grid = np.meshgrid(coordinate, coordinate, coordinate, indexing="ij")
    radius = np.sqrt(sum(axis * axis for axis in grid))
    scale = np.float32(size / 24)
    white = np.where(radius <= 4 * scale, 230.0, 0.0).astype(np.float32)
    gray = np.where((radius > 5 * scale) & (radius <= 9 * scale), 220.0, 0.0).astype(np.float32)
    atlas = np.where((radius > 5 * scale) & (radius <= 9 * scale), 254, 0).astype(np.uint8)
    baseline_rss = _peak_rss_bytes()
    started = time.perf_counter()
    result = ace_enhance_arrays(gray, white, multi_atlas_csf=atlas)
    elapsed = time.perf_counter() - started
    peak_rss = _peak_rss_bytes()
    return {
        "shape": [size, size, size],
        "iterations": 25,
        "synthetic_recipe": "concentric white/gray matter phantom, radii scaled from the committed 24-cubed ACE fixture",
        "elapsed_seconds": round(elapsed, 6),
        "peak_rss_bytes": peak_rss,
        "peak_rss_delta_after_import_and_inputs_bytes": max(0, peak_rss - baseline_rss),
        "enhanced_nonzero_voxels": int(np.count_nonzero(result.enhanced_gray_matter)),
        "cleaned_skeleton_voxels": int(np.count_nonzero(result.skeleton)),
        "thinned_skeleton_voxels": int(np.count_nonzero(result.thinned_skeleton)),
    }


def main() -> None:
    if len(sys.argv) == 3 and sys.argv[1] == "--worker":
        print(json.dumps(_run_size(int(sys.argv[2])), sort_keys=True))
        return
    sizes = (32, 40, 48)
    rows = []
    for size in sizes:
        started = time.perf_counter()
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--worker", str(size)],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        row = json.loads(completed.stdout)
        row["child_wall_seconds"] = round(time.perf_counter() - started, 6)
        rows.append(row)
    print(json.dumps({"benchmark": "synthetic ACE selected V3 defaults", "runs": rows}, indent=2))


if __name__ == "__main__":
    main()
