#!/usr/bin/env python3
"""Freeze the exact iteration-5 checkpoint and registered priors for i6 thinning."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np

from pymacruise.image_io import nifti_ras_on_reference, save_nifti_ras
from pymacruise.lesion_atlas_priors import load_lesion_atlas_shape_priors
from pymacruise.lesion_registered_shapes import transform_registered_lesion_shape_priors
from pymacruise.toads_shapes import ToadShapeBounds


REPO = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = REPO / "build/plus-vandy-atlas-aligned-fine"
DEFAULT_FIXTURE = REPO / "tests/data/synthetic/plus_iteration_6_thinning"
DEFAULT_DESCRIPTOR = REPO / "tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verified(record: dict, label: str) -> Path:
    path = Path(record["path"]).expanduser().resolve()
    if not path.is_file() or path.is_symlink() or sha256(path) != record.get("sha256"):
        raise ValueError(f"{label} is missing or its checkpoint hash does not match: {path}")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--descriptor", type=Path, default=DEFAULT_DESCRIPTOR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIXTURE)
    args = parser.parse_args()
    root = args.workflow_root.resolve()
    descriptor = args.descriptor.resolve()
    destination = args.output_dir.resolve()

    checkpoint_path = root / "loop-chunks/iteration-0005/lesion-registered-loop.json"
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    if checkpoint.get("status") != "complete":
        raise ValueError("iteration-5 registered checkpoint is not complete")
    summary_path = verified(checkpoint["outputs"]["summary"], "iteration-5 summary")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if int(summary.get("completed_iteration", -1)) != 5:
        raise ValueError("expected the saved iteration-5 checkpoint")
    recrop_manifest_path = verified(checkpoint["inputs"]["registered_recrop_manifest"], "registered recrop manifest")
    recrop_manifest = json.loads(recrop_manifest_path.read_text(encoding="utf-8"))
    recrop_summary_path = verified(recrop_manifest["outputs"]["summary"], "registered recrop summary")
    recrop_summary = json.loads(recrop_summary_path.read_text(encoding="utf-8"))

    segmentation_path = verified(checkpoint["outputs"]["segmentation"], "iteration-5 segmentation")
    source_images = {
        name: verified(checkpoint["outputs"][key], f"iteration-5 {name}")
        for name, key in (
            ("skeleton", "skeleton"),
            ("memberships_17ch", "memberships"),
            ("availability_13ch", "availability"),
            ("ordering", "ordering"),
        )
    }
    reference = nib.load(str(segmentation_path))
    if nib.aff2axcodes(reference.affine) != ("R", "A", "S"):
        raise ValueError("the saved iteration-5 segmentation is not in RAS")
    dimensions = reference.shape
    if dimensions != (71, 89, 67):
        raise ValueError(f"unexpected iteration-5 grid: {dimensions}")
    loaded = {"segmentation": np.asarray(reference.dataobj, dtype=np.int16)}
    for name, path in source_images.items():
        image = nib.load(str(path))
        if image.shape[:3] != dimensions or nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise ValueError(f"iteration-5 {name} is not on the checkpoint RAS grid")
        if not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-5):
            raise ValueError(f"iteration-5 {name} has a mismatched affine")
        loaded[name] = np.asarray(image.dataobj)
    if loaded["memberships_17ch"].shape != (*dimensions, 17):
        raise ValueError("expected the 13 class memberships plus four auxiliary channels")
    if loaded["availability_13ch"].shape != (*dimensions, 13):
        raise ValueError("expected 13 class availability channels")

    descriptor_shapes = load_lesion_atlas_shape_priors(descriptor)
    class_names = tuple(recrop_summary["class_names"])
    if class_names != tuple(descriptor_shapes.shape_maps):
        raise ValueError("registered checkpoint and descriptor class order differ")
    prior_maps = tuple(descriptor_shapes.shape_maps[name] for name in class_names)
    bounds = tuple(
        ToadShapeBounds(
            tuple(int(value) for value in row["lower_xyz"]),
            tuple(int(value) for value in row["upper_xyz"]),
            int(row["support_count"]),
            None,
        )
        for row in summary["registration_shape_bounds"]
    )
    if len(bounds) != 13:
        raise ValueError("iteration-5 summary must provide 13 registered shape bounds")
    modality_path = verified(recrop_manifest["outputs"]["modalities"], "registered recrop modalities")
    modality_image = nib.load(str(modality_path))
    if modality_image.shape[:3] != dimensions or not np.allclose(modality_image.affine, reference.affine, rtol=0, atol=1e-5):
        raise ValueError("registered recrop modality geometry does not match iteration-5 checkpoint")
    template_image = nib.load(str(descriptor_shapes.files[class_names[0]]))
    image_resolution = tuple(float(value) for value in nib.affines.voxel_sizes(modality_image.affine))
    template_resolution = tuple(float(value) for value in nib.affines.voxel_sizes(template_image.affine))
    registration = SimpleNamespace(
        transform=np.asarray(summary["registration_transform"], dtype=np.float32),
        shape_bounds=bounds,
        memberships=np.zeros((*dimensions, 13), dtype=np.float32),
    )
    shapes = transform_registered_lesion_shape_priors(
        registration,
        prior_maps,
        image_resolutions_mm=image_resolution,
        template_resolutions_mm=template_resolution,
    )
    saved_shapes_path = verified(recrop_manifest["outputs"]["shapes"], "saved transformed atlas shapes")
    saved_shapes_image = nib.load(str(saved_shapes_path))
    saved_shapes = np.asarray(saved_shapes_image.dataobj, dtype=np.float32)
    if saved_shapes.shape != shapes.shape or not np.array_equal(saved_shapes, shapes):
        raise ValueError("recomputed iteration-5 registered shape priors differ from saved shape volume")
    registered_shape_values_hash = hashlib.sha256(
        np.asarray(shapes, dtype="<f4").ravel(order="C").tobytes()
    ).hexdigest()
    saved_shape_values_hash = hashlib.sha256(
        np.asarray(saved_shapes, dtype="<f4").ravel(order="C").tobytes()
    ).hexdigest()
    if registered_shape_values_hash != saved_shape_values_hash:
        raise ValueError("registered shape-prior voxel bits differ from the saved recrop output")

    arrays = {
        "segmentation": loaded["segmentation"],
        "skeleton": loaded["skeleton"].astype(np.int16, copy=False),
        "memberships_17ch": loaded["memberships_17ch"].astype(np.float32, copy=False),
        "availability_13ch": loaded["availability_13ch"].astype(np.uint8, copy=False),
        "ordering": loaded["ordering"].astype(np.float32, copy=False),
        "registered_shapes_13ch": shapes,
    }
    destination.mkdir(parents=True, exist_ok=True)
    files = {}
    for name, values in arrays.items():
        output = destination / f"{name}.nii.gz"
        save_nifti_ras(nifti_ras_on_reference(reference, values), output, dtype=values.dtype)
        files[name] = {"path": output.name, "sha256": sha256(output), "dtype": str(values.dtype), "shape": list(values.shape)}

    provenance = {
        "case": "P41z synthetic iteration-6 LesionTOADS thinning input",
        "synthetic_only": True,
        "input_checkpoint_iteration": 5,
        "source_order": ["LesionToadSegmentation.initThinningLabels", "LesionToadSegmentation.propagateThinning"],
        "source_call_order_sha256": {
            "LesionToadSegmentation.java": "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99",
            "AlgorithmLesionToads.java": "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb",
        },
        "java_reference_harness": {
            "path": "tests/java_reference/Lesion13ClassNextThinningReferenceHarness.java",
            "sha256": "101aa12bfe874d937b4a0f8f208471dfe38139aff739f2f1d61cda16a159b83f",
            "reference_kind": "selected-method Java excerpt with copied heap and topology LUT helpers",
        },
        "parameters": {
            "classes": 13,
            "lesion_class_index": 10,
            "volume_ratio": 0.5,
            "limit": 1000.0,
            "use_lesion_max": False,
            "class_names": list(class_names),
            "object_types": [item.topology for item in descriptor_shapes.structures],
            "membership_channels": {"classes": [0, 12], "lesion": 13, "black_hole": 14, "outlier": 15, "relation": 16},
            "spatial_convention": "RAS",
            "image_resolutions_mm": list(image_resolution),
            "template_resolutions_mm": list(template_resolution),
            "transformed_shape_comparison": "bit-exact to registered_recrop/transformed_shapes_13.nii.gz",
            "transformed_shape_voxel_values_sha256": registered_shape_values_hash,
            "saved_transformed_shape_voxel_values_sha256": saved_shape_values_hash,
        },
        "source_checkpoint_sha256": {
            "iteration5_manifest": sha256(checkpoint_path),
            "iteration5_summary": sha256(summary_path),
            "registered_recrop_manifest": sha256(recrop_manifest_path),
            "registered_recrop_summary": sha256(recrop_summary_path),
            "atlas_descriptor": sha256(descriptor),
            "atlas_shape_priors": {
                name: sha256(descriptor_shapes.files[name]) for name in class_names
            },
            "iteration5_segmentation": sha256(segmentation_path),
            **{f"iteration5_{name}": sha256(path) for name, path in source_images.items()},
            "registered_recrop_modalities": sha256(modality_path),
            "saved_registered_shapes": sha256(saved_shapes_path),
        },
        "files": files,
    }
    (destination / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"fixture_dir": str(destination), "files": len(files), "transformed_shape_bit_exact": True}, indent=2))


if __name__ == "__main__":
    main()
