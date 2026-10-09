#!/usr/bin/env python3
"""Compare saved plus atlas-boundary output with selected JIST Java source."""
from __future__ import annotations

import argparse
import hashlib
import json
import struct
from pathlib import Path

import numpy as np


def read_int_volume(path: Path, shape: tuple[int, int, int]) -> np.ndarray:
    raw = path.read_bytes()
    found = struct.unpack(">III", raw[:12])
    if found != shape:
        raise SystemExit(f"{path.name}: output shape {found} != {shape}")
    count = int(np.prod(shape))
    values = np.frombuffer(raw, dtype=">i4", offset=12)
    if values.size != count:
        raise SystemExit(f"{path.name}: expected {count} integers, found {values.size}")
    return values.astype(np.int32).reshape(shape)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    root = args.reference_dir
    provenance = json.loads((root / "input-provenance.json").read_text())
    shape = tuple(provenance["shape_xyz"])
    source = Path(__file__).resolve().parents[1] / "tests/java_reference/upstream/edu/vanderbilt/masi/plugins/CRUISE/utilities/MinFilter3D.java"
    import re
    def case(java_name: str, python_name: str, log_name: str, input_hash: str, expected_hash: str) -> dict:
        java = read_int_volume(root / java_name, shape)
        python = np.load(root / python_name).astype(np.int32)
        if python.shape != shape:
            raise SystemExit(f"{python_name}: Python boundary shape {python.shape} != {shape}")
        difference = java != python
        match = re.search(r"peak_sampled_rss_gb=([0-9.]+)", (root / log_name).read_text())
        if match is None:
            raise SystemExit(f"{log_name}: watchdog log lacks peak RSS")
        return {
            "input_sha256": input_hash,
            "java_reference_output_sha256": hashlib.sha256((root / java_name).read_bytes()).hexdigest(),
            "python_reference_output_sha256": expected_hash,
            "java_peak_sampled_rss_gb": float(match.group(1)),
            "voxel_count": int(java.size),
            "java_foreground_voxels": int(np.count_nonzero(java)),
            "python_foreground_voxels": int(np.count_nonzero(python)),
            "changed_voxels": int(np.count_nonzero(difference)),
            "exact": bool(not difference.any()),
        }
    actual = case(
        "java-raw-boundary.bin", "python-raw-boundary.npy", "boundary-java-watchdog.log",
        provenance["boundary_input_binary_sha256"], provenance["raw_boundary_sha256"],
    )
    folded = case(
        "java-folded-raw-boundary.bin", "python-folded-raw-boundary.npy", "boundary-folded-java-watchdog.log",
        provenance["folded_boundary_input_sha256"], provenance["folded_raw_boundary_sha256"],
    )
    report = {
        "synthetic_only": True,
        "comparison": "preserved JIST MinFilter3D source versus Python plus raw atlas boundary on held-fixed saved workflow labels",
        "scope": "The same RAS segmentation crop was passed through selected Java MinFilter3D.GetBoundrayForOneLabel and the Python run_plus_boundary_stage; only the source boundary algorithm is compared.",
        "workflow_run_signature": provenance["workflow_run_signature"],
        "implementation_signature": provenance["implementation_signature"],
        "shape_xyz": list(shape),
        "coordinate_system": "RAS",
        "selected_java_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "java_rss_limit_gb": 8,
        "saved_plus_crop_case": actual,
        "folded_synthetic_case": folded,
        "saved_plus_crop_is_nontrivial": actual["java_foreground_voxels"] > 0,
        "folded_case_is_nontrivial": folded["java_foreground_voxels"] > 0,
        "passed": actual["exact"] and folded["exact"] and folded["java_foreground_voxels"] > 0,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit("plus raw-boundary output differs from selected JIST source or the synthetic case is trivial")


if __name__ == "__main__":
    main()
