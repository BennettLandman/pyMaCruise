#!/usr/bin/env python3
"""Locate per-pass differences in the held-fixed Java/Python TGDM smoother."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from pymacruise.ace import _signed_distance_fmm
from pymacruise.tgdm import _SMOOTHING_PROFILE_SETTINGS, tgdm_evolve


def _read_java(path: Path, shape: tuple[int, int, int]) -> np.ndarray:
    with gzip.open(path, "rb") as stream:
        data = np.frombuffer(stream.read(), dtype=">f4")
    if data.size != int(np.prod(shape)):
        raise SystemExit(f"Java volume has the wrong size: {path}")
    return data.astype(np.float32).reshape(shape)


def _metrics(reference: np.ndarray, candidate: np.ndarray) -> dict[str, float | int]:
    difference = np.abs(reference.astype(np.float64) - candidate.astype(np.float64))
    return {
        "sign_differences": int(np.count_nonzero(np.signbit(reference) != np.signbit(candidate))),
        "changed_voxels": int(np.count_nonzero(reference != candidate)),
        "mean_absolute_voxel_error": float(difference.mean()),
        "maximum_absolute_voxel_error": float(difference.max(initial=0.0)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument(
        "--source-float-rounding", action="store_true",
        help="round fast-marching values to float32 at Java image/heap writes",
    )
    args = parser.parse_args()
    root = args.reference_dir.resolve()
    provenance = json.loads((root / "input-provenance.json").read_text())
    shape = tuple(provenance["shape_xyz"])
    initial = np.load(root / "python-initial_level_set.npy").astype(np.float32)
    smoothing_input = _signed_distance_fmm(
        initial,
        max_distance=8.0,
        round_updates_to_float32=args.source_float_rounding,
    )
    java_smoothing_input = _read_java(root / "java/smoothing_input.raw.gz", shape)
    input_metrics = _metrics(java_smoothing_input, smoothing_input)
    python = tgdm_evolve(
        smoothing_input,
        profile="central",
        curvature_force=-1.0,
        pressure_force=0.0,
        outer_iterations=6,
        settings=_SMOOTHING_PROFILE_SETTINGS,
        round_reinitialization_updates_to_float32=args.source_float_rounding,
        record_passes=True,
    )
    passes = []
    for index, (pre, post) in enumerate(
        zip(python.pre_reinitialization_level_sets, python.pass_level_sets), start=1
    ):
        java_pre = _read_java(root / f"java/smoothed_pre_reinit_pass_{index:02d}.raw.gz", shape)
        java_post = _read_java(root / f"java/smoothed_pass_{index:02d}.raw.gz", shape)
        passes.append({
            "pass": index,
            "before_reinitialization": _metrics(java_pre, pre),
            "after_reinitialization": _metrics(java_post, post),
        })
    first_mismatch = next(
        (
            {"pass": item["pass"], "stage": stage}
            for item in passes
            for stage, key in (
                ("curvature_update", "before_reinitialization"),
                ("fast_marching_reinitialization", "after_reinitialization"),
            )
            if item[key]["sign_differences"] > 0
        ),
        None,
    )
    report = {
        "synthetic_only": True,
        "diagnostic_only": True,
        "comparison": "selected Java SmoothCortexProfile pass outputs vs Python TGDM smoother",
        "scope": "isolates input fast marching, per-pass curvature updates, and per-pass fast-marching reinitialization; does not alter production outputs",
        "workflow_run_signature": provenance["workflow_run_signature"],
        "shape_xyz": list(shape),
        "coordinate_system": "RAS",
        "settings": {
            "outer_iterations": 6,
            "internal_iterations": 30,
            "final_internal_iterations": 30,
            "delta_t": 0.1,
            "final_delta_t": 0.01,
            "max_signed_distance": 6.2,
            "final_max_signed_distance": 20.0,
            "round_distance_updates_to_float32": args.source_float_rounding,
        },
        "initial_reinitialization": input_metrics,
        "passes": passes,
        "first_sign_mismatch": first_mismatch,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
