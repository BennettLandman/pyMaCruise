"""Measure the default TGDM pipeline on bundled synthetic data only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import resource
import sys
import time

import nibabel as nib
import numpy as np

from pymacruise.tgdm import TGDMBaselineParameters, tgdm_nested_arrays


SYNTHETIC = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "tgdm"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--size",
        type=int,
        default=25,
        help="Cubic synthetic phantom size; 25 uses the committed TGDM inputs.",
    )
    args = parser.parse_args()
    if args.size < 17:
        parser.error("--size must be at least 17 voxels")

    if args.size == 25:
        membership = _load("initial_membership")
        wm = _load("white_matter")
        gm = _load("gray_matter")
        ace_gm = _load("ace_gray_matter")
        skeleton = _load("ace_skeleton")
    else:
        membership, wm, gm, ace_gm, skeleton = _make_sphere(args.size)

    start = time.perf_counter()
    result = tgdm_nested_arrays(
        membership,
        wm,
        gm,
        ace_gm,
        ace_skeleton=skeleton,
        gvf=np.zeros((*membership.shape, 3), dtype=np.float32),
        parameters=TGDMBaselineParameters(),
    )
    elapsed = time.perf_counter() - start
    max_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":
        max_rss *= 1024
    report = {
        "synthetic_only": True,
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "shape": list(membership.shape),
        "default_wrapper_parameters": True,
        "synthetic_inputs": "concentric tissues, ACE skeleton, zero GVF, no brain mask",
        "elapsed_seconds": round(elapsed, 3),
        "process_peak_rss_mib": round(max_rss / (1024 * 1024), 1),
        "completed_pde_steps": {
            "initial_smoothing": result.smoothed_initial.completed_steps,
            "inner_seed": result.initial_inner.completed_steps,
            "central": result.central.completed_steps,
            "central_shrink": result.central_shrink.completed_steps,
            "inner_final": result.inner.completed_steps,
            "outer": result.outer.completed_steps,
        },
        "nonpositive_voxel_counts": {
            "inner": int(np.count_nonzero(result.inner.level_set <= 0)),
            "central": int(np.count_nonzero(result.central.level_set <= 0)),
            "outer": int(np.count_nonzero(result.outer.level_set <= 0)),
        },
        "nested_mask_containment": {
            "inner_inside_central": _containment(result.inner.level_set, result.central.level_set),
            "central_inside_outer": _containment(result.central.level_set, result.outer.level_set),
        },
        "output_finite": all(
            np.isfinite(volume).all()
            for volume in (
                result.initial_level_set,
                result.smoothed_initial.level_set,
                result.initial_inner.level_set,
                result.central_shrink.level_set,
                result.inner.level_set,
                result.central.level_set,
                result.outer.level_set,
            )
        ),
    }
    print(json.dumps(report, indent=2))


def _containment(inner: np.ndarray, outer: np.ndarray) -> float:
    inner_mask = np.asarray(inner) <= 0
    outer_mask = np.asarray(outer) <= 0
    count = int(inner_mask.sum())
    return float(np.count_nonzero(inner_mask & outer_mask) / count) if count else 1.0


def _load(name: str) -> np.ndarray:
    return np.asarray(nib.load(str(SYNTHETIC / f"{name}.nii.gz")).dataobj)


def _make_sphere(size: int) -> tuple[np.ndarray, ...]:
    center = (size - 1) / 2
    x, y, z = np.indices((size, size, size), dtype=np.float32)
    radius = np.sqrt((x - center) ** 2 + (y - center) ** 2 + (z - center) ** 2)
    white_radius = size * 0.2
    gray_radius = size * 0.32
    membership = np.where(radius <= gray_radius, 255, 0).astype(np.uint8)
    wm = np.where(radius <= white_radius, 255, 0).astype(np.uint8)
    gray = np.where((radius > white_radius) & (radius <= gray_radius), 220, 0).astype(np.uint8)
    ace_gm = gray.copy()
    skeleton = np.zeros_like(membership)
    skeleton[
        (np.abs(radius - white_radius - 0.75) < 0.5) & (z == int(center))
    ] = 254
    return membership, wm, gray, ace_gm, skeleton


if __name__ == "__main__":
    main()
