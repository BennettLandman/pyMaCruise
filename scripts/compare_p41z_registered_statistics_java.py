#!/usr/bin/env python3
"""Compare selected Java statistics with a saved P41z registered state.

This checks the VandyCRUISE computeCentroids/computeVariances/
computeChannelWeights excerpt on the exact registered state produced by a
persisted synthetic LesionTOADS loop. It does not run thinning or compare a
complete registered-loop iteration.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from compare_p41z_next_iteration_statistics_java import compare  # noqa: E402


def _checked_path(record: dict, *, name: str) -> Path:
    path = Path(record.get("path", "")).resolve()
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(f"registered state is missing {name}: {path}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != record.get("sha256"):
        raise ValueError(f"registered state {name} fails its checkpoint hash")
    return path


def _load_record(manifest: dict, name: str) -> tuple[Path, nib.spatialimages.SpatialImage]:
    record = manifest["outputs"][name]
    path = _checked_path(record, name=name)
    image = nib.load(str(path))
    if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
        raise ValueError(f"registered state {name} must be RAS")
    return path, image


def prepare_case(loop_manifest_path: Path, output_dir: Path) -> Path:
    loop_path = loop_manifest_path.expanduser().resolve()
    loop = json.loads(loop_path.read_text(encoding="utf-8"))
    if loop.get("status") != "complete" or loop.get("stage") != "plus-lesion-registered-loop":
        raise ValueError("input must be a completed plus registered-loop checkpoint")
    if loop.get("spatial_convention") != "RAS":
        raise ValueError("registered-loop checkpoint must be RAS")
    recrop_record = loop["inputs"]["registered_recrop_manifest"]
    recrop_path = _checked_path(recrop_record, name="registered recrop manifest")
    recrop = json.loads(recrop_path.read_text(encoding="utf-8"))
    image_record = recrop["outputs"]["modalities"]
    image_path = _checked_path(image_record, name="registered modalities")
    images = nib.load(str(image_path))
    if len(images.shape) != 4 or images.shape[-1] != 2 or nib.aff2axcodes(images.affine) != ("R", "A", "S"):
        raise ValueError("registered modalities must be two-channel RAS NIfTI")

    membership_path, membership_image = _load_record(loop, "memberships")
    segmentation_path, segmentation_image = _load_record(loop, "segmentation")
    if membership_image.shape != (*images.shape[:3], 17) or segmentation_image.shape != images.shape[:3]:
        raise ValueError("registered modalities, memberships, and segmentation grids differ")
    if not np.allclose(images.affine, membership_image.affine, rtol=0, atol=1e-4) or not np.allclose(
        images.affine, segmentation_image.affine, rtol=0, atol=1e-4
    ):
        raise ValueError("registered modalities, memberships, and segmentation affines differ")

    summary_path = _checked_path(loop["outputs"]["summary"], name="registered summary")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    centroids = summary["refreshed_centroids"]
    initial_case = json.loads(
        (ROOT / "tests/data/synthetic/p41z_iteration6_statistics_java/case.json").read_text(
            encoding="utf-8"
        )
    )
    case = {
        "case": "P41z registered-state selected Java statistics comparison",
        "modalities": initial_case["modalities"],
        "object_types": initial_case["object_types"],
        "initial_class_centroids": centroids["class_centroids"],
        "initial_lesion_centroids": centroids["lesion_centroids"],
        "initial_black_hole_centroids": centroids["black_hole_centroids"],
        "tissue_power_exponents": initial_case["tissue_power_exponents"],
        "lesion_power_exponents": initial_case["lesion_power_exponents"],
        "norm_type": initial_case["norm_type"],
        "outlier_distance": initial_case["outlier_distance"],
        "optimized_csf_all_modalities": initial_case["optimized_csf_all_modalities"],
        "input_files": {
            "images": "images_t1_flair.nii.gz",
            "memberships": "class_memberships_13ch.nii.gz",
            "segmentation": "registered_segmentation.nii.gz",
            "lesion_membership": "lesion_membership.nii.gz",
            "black_hole_membership": "black_hole_membership.nii.gz",
            "outlier_membership": "outlier_membership.nii.gz",
        },
        "provenance": {
            "registered_loop_manifest": str(loop_path),
            "registered_recrop_manifest": str(recrop_path),
            "registered_modalities": str(image_path),
            "registered_memberships": str(membership_path),
            "registered_segmentation": str(segmentation_path),
            "completed_iteration": summary.get("completed_iteration"),
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    affine = images.affine
    arrays = {
        "images_t1_flair.nii.gz": np.asarray(images.dataobj, dtype=np.float32),
        "class_memberships_13ch.nii.gz": np.asarray(membership_image.dataobj[..., :13], dtype=np.float32),
        "registered_segmentation.nii.gz": np.asarray(segmentation_image.dataobj, dtype=np.int16),
        "lesion_membership.nii.gz": np.asarray(membership_image.dataobj[..., 13], dtype=np.float32),
        "black_hole_membership.nii.gz": np.asarray(membership_image.dataobj[..., 14], dtype=np.float32),
        "outlier_membership.nii.gz": np.asarray(membership_image.dataobj[..., 15], dtype=np.float32),
    }
    for name, values in arrays.items():
        nib.save(nib.Nifti1Image(values, affine), str(output_dir / name))
    case_path = output_dir / "case.json"
    case_path.write_text(json.dumps(case, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return case_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("loop_manifest", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    case_path = prepare_case(args.loop_manifest, args.output_dir)
    result = compare(case_path.parent)
    result["case_provenance"] = json.loads(case_path.read_text(encoding="utf-8"))["provenance"]
    report = args.report.resolve()
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
