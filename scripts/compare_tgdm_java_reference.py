"""Compare the default Python TGDM pipeline with committed Java phantoms."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import sys

import nibabel as nib
import numpy as np
from skimage.measure import marching_cubes

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from pymacruise.tgdm import tgdm_nested_images  # noqa: E402
from pymacruise.validation import (  # noqa: E402
    dice_coefficient,
    mesh_surface_distance_metrics,
    voxel_error_metrics,
    write_metrics_report,
)

DATA_ROOT = REPO_ROOT / "tests" / "data" / "synthetic" / "tgdm"
REFERENCE_ROOT = DATA_ROOT / "java_reference"
LEVEL_SET_FIELDS = {
    "signed_seed": "initial_level_set",
    "smoothed_initial": "smoothed_initial",
    "initial_inner": "initial_inner",
    "central": "central",
    "central_shrink": "central_shrink",
    "inner": "inner",
    "outer": "outer",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, help="also write this JSON report")
    args = parser.parse_args()

    inputs = {
        name: nib.load(str(DATA_ROOT / f"{name}.nii.gz"))
        for name in (
            "initial_membership",
            "white_matter",
            "gray_matter",
            "ace_gray_matter",
            "ace_skeleton",
        )
    }
    gvf = nib.Nifti1Image(
        np.zeros((*inputs["initial_membership"].shape, 3), dtype=np.float32),
        inputs["initial_membership"].affine,
    )
    result = tgdm_nested_images(
        inputs["initial_membership"],
        inputs["white_matter"],
        inputs["gray_matter"],
        inputs["ace_gray_matter"],
        ace_skeleton=inputs["ace_skeleton"],
        gvf=gvf,
    )

    stages: dict[str, object] = {}
    for reference_name, result_name in LEVEL_SET_FIELDS.items():
        reference_image = nib.load(str(REFERENCE_ROOT / f"{reference_name}.nii.gz"))
        result_value = getattr(result.execution, result_name)
        candidate_volume = (
            result_value if result_name == "initial_level_set" else result_value.level_set
        )
        candidate = np.asarray(candidate_volume, dtype=np.float32)
        reference = np.asarray(reference_image.dataobj, dtype=np.float32)
        assert_ras_identity(reference_image)
        candidate_vertices, candidate_faces, _, _ = marching_cubes(candidate, level=0.0)
        reference_vertices, reference_faces, _, _ = marching_cubes(reference, level=0.0)
        stages[reference_name] = {
            "dice": dice_coefficient(candidate <= 0, reference <= 0),
            **voxel_error_metrics(candidate, reference),
            **mesh_surface_distance_metrics(
                candidate_vertices,
                candidate_faces,
                reference_vertices,
                reference_faces,
            ),
        }

    tissue_reference = nib.load(str(REFERENCE_ROOT / "tissue_sign.nii.gz"))
    assert_ras_identity(tissue_reference)
    tissue_candidate = np.asarray(result.tissue_sign.dataobj)
    tissue_expected = np.asarray(tissue_reference.dataobj)
    stages["tissue_sign"] = {
        "exact": bool(np.array_equal(tissue_candidate, tissue_expected)),
        **voxel_error_metrics(tissue_candidate, tissue_expected),
    }

    report = {
        "synthetic_only": True,
        "reference": "MaGenericTGDM CVS 1.5; DistanceFieldyh CVS 1.1",
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "nibabel_version": nib.__version__,
        "orientation": "RAS",
        "voxel_spacing_mm": [1.0, 1.0, 1.0],
        "acceptance": {
            "per_stage_mask_dice_minimum": 0.97,
            "symmetric_mean_surface_distance_mm_maximum": 0.01,
            "symmetric_hausdorff_mm_maximum": 1.0,
        },
        "stages": stages,
    }
    report["passed"] = all(
        metrics["dice"] >= 0.97
        and metrics["mean_absolute_distance_mm"] <= 0.01
        and metrics["hausdorff_mm"] <= 1.0
        for stage, metrics in stages.items()
        if stage != "tissue_sign"
    ) and stages["tissue_sign"]["exact"]
    serialized = json.dumps(report, indent=2, sort_keys=True)
    print(serialized)
    if args.report:
        write_metrics_report(args.report, report)
    if not report["passed"]:
        raise SystemExit("TGDM synthetic reference comparison did not meet acceptance thresholds")


def assert_ras_identity(image: nib.spatialimages.SpatialImage) -> None:
    if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
        raise ValueError("Java TGDM reference output must be RAS-oriented")
    if not np.allclose(image.affine, np.eye(4), atol=0.0):
        raise ValueError("TGDM synthetic reference must use its identity RAS affine")


if __name__ == "__main__":
    main()
