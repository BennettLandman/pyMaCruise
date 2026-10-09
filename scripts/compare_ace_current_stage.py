"""Compare selected Java ACE outputs with a saved Python CRUISE-core stage."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


OUTPUTS = {
    "enhanced_gray_matter": ("enhanced_gm_java.nii.gz", "ace_enhanced_gm.nii.gz"),
    "skeleton": ("csf_skeleton_java.nii.gz", "ace_skeleton.nii.gz"),
    "thinned_skeleton": ("thinned_skeleton_java.nii.gz", "ace_thinned_skeleton.nii.gz"),
}


def compare(reference_dir: Path, stage_dir: Path, diffusion_iterations: int) -> dict[str, object]:
    outputs: dict[str, object] = {}
    for name, (java_filename, python_filename) in OUTPUTS.items():
        java_path = reference_dir / java_filename
        python_path = stage_dir / python_filename
        java_image = nib.load(str(java_path))
        python_image = nib.load(str(python_path))
        if nib.aff2axcodes(java_image.affine) != ("R", "A", "S"):
            raise ValueError(f"Java reference is not RAS: {java_filename}")
        if nib.aff2axcodes(python_image.affine) != ("R", "A", "S"):
            raise ValueError(f"Python workflow result is not RAS: {python_filename}")
        if java_image.shape != python_image.shape or not np.allclose(
            java_image.affine, python_image.affine, atol=1e-5
        ):
            raise ValueError(f"Java and Python outputs do not share one RAS grid: {name}")
        java = np.asarray(java_image.dataobj)
        python = np.asarray(python_image.dataobj)
        if not np.isfinite(java).all() or not np.isfinite(python).all():
            raise ValueError(f"non-finite ACE output: {name}")
        delta = np.abs(java.astype(np.float64) - python.astype(np.float64))
        record: dict[str, object] = {
            "exact": bool(np.array_equal(java, python)),
            "changed_voxels": int(np.count_nonzero(java != python)),
            "max_absolute_error": float(delta.max(initial=0.0)),
            "mean_absolute_error": float(delta.mean()),
            "java_sha256": hashlib.sha256(java_path.read_bytes()).hexdigest(),
            "python_sha256": hashlib.sha256(python_path.read_bytes()).hexdigest(),
        }
        if "skeleton" in name:
            java_mask = java > 0
            python_mask = python > 0
            intersection = int(np.count_nonzero(java_mask & python_mask))
            denominator = int(np.count_nonzero(java_mask) + np.count_nonzero(python_mask))
            record.update(
                {
                    "java_foreground": int(np.count_nonzero(java_mask)),
                    "python_foreground": int(np.count_nonzero(python_mask)),
                    "dice": 2.0 * intersection / denominator if denominator else 1.0,
                }
            )
        outputs[name] = record
    return {
        "synthetic_only": True,
        "orientation": "RAS",
        "shape": list(nib.load(str(stage_dir / "ace_enhanced_gm.nii.gz")).shape),
        "diffusion_iterations": diffusion_iterations,
        "comparison": "selected Java ACE output vs current Python CRUISE-core outputs",
        "outputs": outputs,
        "all_outputs_exact": all(bool(row["exact"]) for row in outputs.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("java_reference_dir", type=Path)
    parser.add_argument("python_stage_dir", type=Path)
    parser.add_argument("--diffusion-iterations", type=int, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.diffusion_iterations < 0:
        parser.error("--diffusion-iterations must be non-negative")
    report = compare(args.java_reference_dir, args.python_stage_dir, args.diffusion_iterations)
    serialized = json.dumps(report, indent=2, sort_keys=True)
    print(serialized)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(serialized + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
