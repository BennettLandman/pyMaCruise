#!/usr/bin/env python3
"""Compare saved human source TOADS memberships with a PyMaCRUISE run.

This produces index-aligned diagnostics only. The source TOADS maps have
zero-origin LPS headers, so this script does not claim physical-grid parity.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import apply_orientation, io_orientation, ornt_transform


SOURCE_OUTPUT = Path("/Users/landmaba/dev/mipav/sample-macruise/output")
PYTHON_RUN = Path("build/invivo-base-orientation-only-source-profile-2026-10-12")
BET_MASK = Path("/Users/landmaba/dev/mipav/sample-macruise/input/T1_bet_mask.nii.gz")
MASEG = Path("/Users/landmaba/dev/mipav/sample-macruise/input/orig_target_seg.nii.gz")


def image(path: Path) -> tuple[nib.spatialimages.SpatialImage, np.ndarray]:
    nii = nib.load(str(path))
    return nii, np.asarray(nii.dataobj)


def dice(a: np.ndarray, b: np.ndarray) -> float:
    a, b = a.astype(bool), b.astype(bool)
    denom = int(a.sum() + b.sum())
    return 1.0 if denom == 0 else float(2 * np.count_nonzero(a & b) / denom)


def mask_record(mask_path: Path, maseg_path: Path, source_mask_path: Path) -> dict[str, object]:
    mask_nii, mask = image(mask_path)
    seg_nii, seg = image(maseg_path)
    source_nii, source = image(source_mask_path)
    for name, candidate in (("MaSeg", seg_nii), ("source MAstrip mask", source_nii)):
        if candidate.shape != mask_nii.shape or not np.allclose(
            candidate.affine, mask_nii.affine, rtol=0, atol=1e-4
        ):
            raise ValueError(f"{name} does not share the supplied brain-mask grid")
    return {
        "geometry_scope": "same raw-T1 grid",
        "supplied_mask_voxels": int(np.count_nonzero(mask)),
        "source_maseg_mask_voxels": int(np.count_nonzero(seg)),
        "source_saved_mask_voxels": int(np.count_nonzero(source)),
        "maseg_binarization_matches_saved_source_mask": bool(
            np.array_equal(seg > 0, source > 0)
        ),
        "supplied_brain_mask_vs_source_maseg_mask_dice": dice(mask > 0, seg > 0),
        "different_voxels": int(np.count_nonzero((mask > 0) != (seg > 0))),
    }


def membership_record(source_path: Path, python_path: Path, threshold: float) -> dict[str, object]:
    source_nii, source = image(source_path)
    python_nii, python = image(python_path)
    transform = ornt_transform(
        io_orientation(source_nii.affine), io_orientation(python_nii.affine)
    )
    source_aligned = apply_orientation(source, transform)
    if source_aligned.shape != python.shape:
        raise ValueError(
            f"shape mismatch after axis orientation: {source_aligned.shape} vs {python.shape}"
        )
    source_f = source_aligned.astype(np.float64, copy=False)
    python_f = python.astype(np.float64, copy=False)
    sample_source, sample_python = source_f.ravel()[::100], python_f.ravel()[::100]
    correlation = float(np.corrcoef(sample_source, sample_python)[0, 1])
    source_active, python_active = source_f > threshold, python_f > threshold
    return {
        "comparison_scope": "voxel indices after axis orientation only",
        "source_axis_codes": list(nib.aff2axcodes(source_nii.affine)),
        "python_axis_codes": list(nib.aff2axcodes(python_nii.affine)),
        "aligned_shape": list(source_aligned.shape),
        "mean_absolute_difference": float(np.mean(np.abs(source_f - python_f))),
        "sampled_correlation_stride_100": correlation,
        "threshold": threshold,
        "threshold_dice": dice(source_active, python_active),
        "source_voxels_above_threshold": int(source_active.sum()),
        "python_voxels_above_threshold": int(python_active.sum()),
        "source_affine": source_nii.affine.tolist(),
        "python_affine": python_nii.affine.tolist(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-output", type=Path, default=SOURCE_OUTPUT)
    parser.add_argument("--python-run", type=Path, default=PYTHON_RUN)
    parser.add_argument("--brain-mask", type=Path, default=BET_MASK)
    parser.add_argument("--maseg", type=Path, default=MASEG)
    parser.add_argument("--threshold", type=float, default=102.0)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    names = ("gm", "wmfill", "wmmask")
    report: dict[str, object] = {
        "scope": "index-aligned human output diagnostic; not physical-space equivalence",
        "source_python_comparison": {},
        "brain_mask_comparison": mask_record(
            args.brain_mask,
            args.maseg,
            args.source_output / "ToadsSeperate/MultiAtlasStripping/T1_N4Corrected_MAstrip_mask.nii.gz",
        ),
    }
    comparisons = report["source_python_comparison"]
    assert isinstance(comparisons, dict)
    for pass_name, source_prefix, python_stage, python_prefix in (
        ("first", "Toads1stPass/T1_N4Corrected_MAstripped_toads1", "toads-first", "t1-n4"),
        ("second", "Toads2ndPass/T1_N4Corrected_MAstripped_durastripped_toads2", "toads-second", "stripped_t1"),
    ):
        comparisons[pass_name] = {}
        pass_comparisons = comparisons[pass_name]
        assert isinstance(pass_comparisons, dict)
        for product in names:
            source_path = args.source_output / f"ToadsSeperate/{source_prefix}{product}.nii.gz"
            python_path = (
                args.python_run
                / f"stages/{python_stage}/attempt-01/dura-inputs/{python_prefix}_{product}.nii.gz"
            )
            pass_comparisons[product] = membership_record(
                source_path, python_path, args.threshold
            )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "report": str(args.report),
        "brain_mask": report["brain_mask_comparison"],
        "source_python_dice": {
            pass_name: {
                name: round(record["threshold_dice"], 6)
                for name, record in pass_data.items()
            }
            for pass_name, pass_data in comparisons.items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
