#!/usr/bin/env python3
"""Compare supplied historical BASE intermediates with a PyMaCRUISE run.

Source outputs are reoriented to closest RAS using only axis permutations and
flips. No interpolation or oblique-affine resampling is performed. The source
tree has no run manifest, so these comparisons are diagnostic and are not
treated as controlled Java/MATLAB source equivalence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compare_pair(name: str, source_path: Path, python_path: Path, kind: str) -> dict[str, object]:
    source_image = nib.load(str(source_path))
    python_image = nib.load(str(python_path))
    if source_image.shape != python_image.shape:
        raise ValueError(
            f"{name}: source shape {source_image.shape} differs from Python {python_image.shape}"
        )
    source_ras = nib.as_closest_canonical(source_image)
    if nib.aff2axcodes(source_ras.affine) != ("R", "A", "S"):
        raise ValueError(f"{name}: source cannot be permuted/flipped to RAS")
    if nib.aff2axcodes(python_image.affine) != ("R", "A", "S"):
        raise ValueError(f"{name}: Python output is not RAS")
    source = np.asarray(source_ras.dataobj)
    python = np.asarray(python_image.dataobj)
    if kind == "float":
        a = source.astype(np.float32, copy=False)
        b = python.astype(np.float32, copy=False)
        delta = np.abs(a.astype(np.float64) - b.astype(np.float64))
        metrics: dict[str, object] = {
            "exact": bool(np.array_equal(a, b)),
            "mae": float(delta.mean()),
            "rmse": float(np.sqrt(np.mean(np.square(delta)))),
            "max_abs": float(delta.max(initial=0.0)),
            "source_range": [float(a.min()), float(a.max())],
            "python_range": [float(b.min()), float(b.max())],
        }
    elif kind == "mask":
        a = source > 0
        b = python > 0
        denominator = int(a.sum()) + int(b.sum())
        metrics = {
            "dice": 1.0 if denominator == 0 else 2.0 * int(np.count_nonzero(a & b)) / denominator,
            "differing_voxels": int(np.count_nonzero(a != b)),
            "source_positive_voxels": int(a.sum()),
            "python_positive_voxels": int(b.sum()),
        }
    elif kind == "labels":
        a = source.astype(np.int64, copy=False)
        b = python.astype(np.int64, copy=False)
        labels = sorted(set(np.unique(a).tolist()) | set(np.unique(b).tolist()))
        per_label: dict[str, float] = {}
        for label in labels:
            if label == 0:
                continue
            am = a == label
            bm = b == label
            denominator = int(am.sum()) + int(bm.sum())
            per_label[str(label)] = (
                1.0 if denominator == 0 else 2.0 * int(np.count_nonzero(am & bm)) / denominator
            )
        scores = list(per_label.values())
        metrics = {
            "differing_voxels": int(np.count_nonzero(a != b)),
            "voxel_count": int(a.size),
            "mean_nonzero_label_dice": float(np.mean(scores)) if scores else 1.0,
            "minimum_nonzero_label_dice": float(min(scores)) if scores else 1.0,
            "per_label_dice": per_label,
        }
    else:
        raise ValueError(f"Unknown comparison kind: {kind}")
    return {
        "kind": kind,
        "source": {
            "path": str(source_path.resolve()),
            "sha256": sha256(source_path),
            "axis_codes_before": list(nib.aff2axcodes(source_image.affine)),
            "axis_codes_after": list(nib.aff2axcodes(source_ras.affine)),
            "affine_after_axis_reorientation": source_ras.affine.tolist(),
        },
        "python": {
            "path": str(python_path.resolve()),
            "sha256": sha256(python_path),
            "axis_codes": list(nib.aff2axcodes(python_image.affine)),
            "affine": python_image.affine.tolist(),
        },
        "shape_ras": list(source_ras.shape),
        "reorientation": "source axis permutations/flips to closest RAS; no interpolation",
        "metrics": metrics,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_output", type=Path)
    parser.add_argument("python_run", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    source = args.source_output.resolve()
    run = args.python_run.resolve()
    pairs = [
        ("n4", source / "ToadsSeperate/N4Correct/T1_N4Corrected.nii.gz",
         run / "stages/n4/attempt-01/t1-n4.nii.gz", "float"),
        ("mask-selection", source / "ToadsSeperate/MultiAtlasStripping/T1_N4Corrected_MAstrip_mask.nii.gz",
         run / "resources/orientation-only-inputs/BrainMask.nii.gz", "mask"),
        ("toads-first-classification", source / "ToadsSeperate/Toads1stPass/T1_N4Corrected_MAstripped_toads1seg.nii.gz",
         run / "stages/toads-first/attempt-01/classification.nii.gz", "labels"),
        ("toads-first-wm-mask", source / "ToadsSeperate/Toads1stPass/T1_N4Corrected_MAstripped_toads1wmmask.nii.gz",
         run / "stages/toads-first/attempt-01/dura-inputs/t1-n4_wmmask.nii.gz", "mask"),
        ("dura-stripped-t1", source / "ToadsSeperate/removeDura/T1_N4Corrected_MAstripped_durastripped.nii.gz",
         run / "stages/dura/attempt-01/stripped_t1.nii.gz", "float"),
        ("dura-brain-mask", source / "ToadsSeperate/removeDura/T1_N4Corrected_MAstripped_duramask.nii.gz",
         run / "stages/dura/attempt-01/brain_mask.nii.gz", "mask"),
        ("toads-second-classification", source / "ToadsSeperate/Toads2ndPass/T1_N4Corrected_MAstripped_durastripped_toads2seg.nii.gz",
         run / "stages/toads-second/attempt-01/classification.nii.gz", "labels"),
        ("toads-second-wm-mask", source / "ToadsSeperate/Toads2ndPass/T1_N4Corrected_MAstripped_durastripped_toads2wmmask.nii.gz",
         run / "stages/toads-second/attempt-01/dura-inputs/stripped_t1_wmmask.nii.gz", "mask"),
    ]
    result: dict[str, object] = {
        "evidence_category": "diagnostic comparison to supplied historical output; not controlled source equivalence",
        "source_run_provenance": "no source run log or manifest was found",
        "reorientation": "permute/flip historical source images to RAS; discard oblique direction cosines; no interpolation",
        "comparisons": {},
    }
    comparisons = result["comparisons"]
    assert isinstance(comparisons, dict)
    for name, source_path, python_path, kind in pairs:
        comparisons[name] = compare_pair(name, source_path, python_path, kind)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {args.report}; {len(comparisons)} RAS-axis-only comparisons")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
