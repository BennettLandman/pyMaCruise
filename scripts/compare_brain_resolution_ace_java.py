"""Compare full-crop Python MaACE outputs with committed selected-Java goldens."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.ace import ace_enhance_arrays


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/data/synthetic/base_workflow_brain_resolution/ace_maace"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    gray_image = nib.load(FIXTURE / "inputs/gray_matter.nii.gz")
    white_image = nib.load(FIXTURE / "inputs/white_matter.nii.gz")
    csf_image = nib.load(FIXTURE / "inputs/multi_atlas_csf.nii.gz")
    for label, image in (("gray matter", gray_image), ("white matter", white_image), ("multi-atlas CSF", csf_image)):
        if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise ValueError(f"{label} input is not RAS-oriented")
        if image.shape != gray_image.shape or not np.allclose(image.affine, gray_image.affine):
            raise ValueError(f"{label} input does not share the gray-matter grid")

    result = ace_enhance_arrays(
        np.asarray(gray_image.dataobj, dtype=np.float32),
        np.asarray(white_image.dataobj, dtype=np.float32),
        multi_atlas_csf=np.asarray(csf_image.dataobj, dtype=np.float32),
        intensity_threshold=126.9,
        tensor_alpha=1.6,
        gaussian_alpha=1.8,
        diffusion_iterations=25,
        use_updated_csf_in_ace=0,
        merge_method="Merge",
    )
    candidates = {
        "enhanced_gm_java.nii.gz": result.enhanced_gray_matter,
        "csf_skeleton_java.nii.gz": result.skeleton,
        "thinned_skeleton_java.nii.gz": result.thinned_skeleton,
    }
    comparisons = {}
    for name, candidate in candidates.items():
        reference_path = FIXTURE / "java_reference" / name
        reference_image = nib.load(reference_path)
        reference = np.asarray(reference_image.dataobj)
        candidate = np.asarray(candidate)
        if candidate.shape != reference.shape:
            raise ValueError(f"{name}: Python shape {candidate.shape} != Java shape {reference.shape}")
        difference = np.abs(candidate.astype(np.float64) - reference.astype(np.float64))
        equal = bool(np.array_equal(candidate, reference))
        comparisons[name] = {
            "exact": equal,
            "changed_voxels": int(np.count_nonzero(candidate != reference)),
            "max_absolute_error": float(difference.max(initial=0)),
            "mean_absolute_error": float(difference.mean()),
            "foreground_java": int(np.count_nonzero(reference)) if "skeleton" in name else None,
            "foreground_python": int(np.count_nonzero(candidate)) if "skeleton" in name else None,
            "java_sha256": sha256(reference_path),
        }
    report = {
        "synthetic_only": True,
        "grid_shape": list(gray_image.shape),
        "orientation": "RAS",
        "parameters": {
            "intensity_threshold": 126.9,
            "tensor_alpha": 1.6,
            "gaussian_alpha": 1.8,
            "diffusion_iterations": 25,
            "place_to_add_csf": 0,
        },
        "outputs": comparisons,
        "passed": all(item["exact"] for item in comparisons.values()),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["passed"]:
        raise SystemExit("brain-resolution MaACE Java comparison differs")


if __name__ == "__main__":
    main()
