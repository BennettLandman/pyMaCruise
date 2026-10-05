"""File-based selected LesionTOADS initial atlas-alignment checkpoint."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras, nifti_ras_on_reference, save_nifti_ras
from pymacruise.lesion_atlas_priors import load_lesion_atlas_topology
from pymacruise.lesion_initial_atlas_alignment import initialize_lesion_atlas_alignment_from_descriptor


@dataclass(frozen=True)
class PlusLesionAtlasInitResult:
    """Persisted membership initialization and initial rigid-atlas transform."""

    manifest_path: Path
    memberships_path: Path
    summary_path: Path
    transform: tuple[float, float, float, float, float, float]
    resumed: bool


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _signature(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _write_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _spacing(image: nib.spatialimages.SpatialImage) -> tuple[float, float, float]:
    values = np.linalg.norm(np.asarray(image.affine[:3, :3], dtype=np.float64), axis=0)
    if not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError("image geometry must have finite positive RAS voxel spacing")
    return tuple(float(value) for value in values)


def _valid_cached(
    manifest: dict[str, Any], signature: str, root: Path, input_image: nib.spatialimages.SpatialImage
) -> PlusLesionAtlasInitResult | None:
    if manifest.get("status") != "complete" or manifest.get("run_signature") != signature:
        return None
    records = manifest.get("outputs", {})
    memberships_record = records.get("memberships")
    summary_record = records.get("summary")
    if not isinstance(memberships_record, dict) or not isinstance(summary_record, dict):
        return None
    memberships_path = Path(memberships_record.get("path", ""))
    summary_path = Path(summary_record.get("path", ""))
    for path, record in ((memberships_path, memberships_record), (summary_path, summary_record)):
        if not path.is_relative_to(root) or not path.is_file() or path.is_symlink():
            return None
        if _sha256(path) != record.get("sha256"):
            return None
    image = nib.load(str(memberships_path))
    if (
        not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image))
        or len(image.shape) != 4
        or image.shape[:3] != input_image.shape[:3]
        or image.shape[-1] != 13
        or nib.aff2axcodes(image.affine) != ("R", "A", "S")
        or not np.allclose(image.affine, input_image.affine, rtol=0, atol=1e-4)
        or not np.isfinite(np.asarray(image.dataobj)).all()
    ):
        return None
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        transform = tuple(float(value) for value in summary["registration_transform"])
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None
    if len(transform) != 6 or not np.isfinite(transform).all():
        return None
    return PlusLesionAtlasInitResult(
        Path(manifest["manifest_path"]), memberships_path, summary_path, transform, True
    )


def run_plus_lesion_initial_atlas_alignment_checkpoint(
    cropped_modalities_path: str | Path,
    crop_manifest_path: str | Path,
    atlas_descriptor_path: str | Path,
    output_dir: str | Path,
    *,
    init_align_iterations: int = 10,
    pyramid_levels: int = 2,
    resume: bool = False,
) -> PlusLesionAtlasInitResult:
    """Initialize 13-class T1/FLAIR memberships and rigid atlas alignment.

    Inputs come from ``plus-crop-lesion-inputs`` and its manifest. Crop bounds
    are stored inclusive there; the source-initial-membership adapter consumes
    half-open bounds, so this handoff advances each high coordinate by one.
    The wrapper's initialization uses the crop-adjusted preprocessing
    transform, descriptor-order intensity/shape priors, T1 then FLAIR channel
    order, and the selected rigid pyramid registration. Output memberships are
    channel-last 4-D RAS NIfTI; the transform and source call order are kept in
    a JSON summary.
    """
    if isinstance(init_align_iterations, bool) or not isinstance(init_align_iterations, int) or init_align_iterations < 0:
        raise ValueError("init_align_iterations must be a nonnegative integer")
    if isinstance(pyramid_levels, bool) or not isinstance(pyramid_levels, int) or pyramid_levels < 1:
        raise ValueError("pyramid_levels must be a positive integer")
    if not isinstance(resume, bool):
        raise TypeError("resume must be boolean")
    cropped_path = Path(cropped_modalities_path).expanduser().resolve()
    crop_manifest = Path(crop_manifest_path).expanduser().resolve()
    descriptor = Path(atlas_descriptor_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    for role, path in (("cropped modalities", cropped_path), ("crop manifest", crop_manifest), ("atlas descriptor", descriptor)):
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Missing regular {role}: {path}")
    if not cropped_path.name.lower().endswith(".nii.gz"):
        raise ValueError("cropped modalities must be .nii.gz")

    crop_data = json.loads(crop_manifest.read_text(encoding="utf-8"))
    if crop_data.get("status") != "complete":
        raise ValueError("crop manifest must describe a complete checkpoint")
    crop_record = crop_data.get("outputs", {}).get("cropped_modalities", {})
    if (
        Path(crop_record.get("path", "")).resolve() != cropped_path
        or _sha256(cropped_path) != crop_record.get("sha256")
    ):
        raise ValueError("cropped modalities do not match the crop manifest output")
    try:
        inclusive_bounds = tuple(int(value) for value in crop_data["crop"]["bounds"])
        padding = int(crop_data["settings"]["padding"])
        transform = np.asarray(crop_data["crop"]["post_crop_transform"], dtype=np.float32)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("crop manifest is missing bounds, padding, or post-crop transform") from exc
    if len(inclusive_bounds) != 6 or transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("crop manifest contains invalid LesionTOADS geometry")
    half_open_bounds = (
        (inclusive_bounds[0], inclusive_bounds[1] + 1),
        (inclusive_bounds[2], inclusive_bounds[3] + 1),
        (inclusive_bounds[4], inclusive_bounds[5] + 1),
    )

    images = load_nifti_ras(cropped_path)
    if len(images.shape) != 4 or images.shape[-1] != 2:
        raise ValueError("cropped modalities must be 4-D channel-last T1 then FLAIR")
    if not np.isfinite(np.asarray(images.dataobj)).all():
        raise ValueError("cropped modalities contain nonfinite values")
    channels_first = np.asarray(images.dataobj, dtype=np.float32).transpose(3, 0, 1, 2)
    if channels_first.shape[1:] != tuple(
        inclusive_bounds[2 * axis + 1] - inclusive_bounds[2 * axis] + 1 + 2 * padding
        for axis in range(3)
    ):
        raise ValueError("cropped modality dimensions do not match the crop manifest")

    atlas_topology = load_lesion_atlas_topology(descriptor)
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "lesion-initial-atlas-alignment.json"
    signature = _signature(
        {
            "stage": "plus-lesion-initial-atlas-alignment",
            "inputs": {
                "cropped_modalities": {"path": str(cropped_path), "sha256": _sha256(cropped_path)},
                "crop_manifest": {"path": str(crop_manifest), "sha256": _sha256(crop_manifest)},
                "atlas_descriptor": {"path": str(descriptor), "sha256": _sha256(descriptor)},
            },
            "settings": {
                "init_align_iterations": init_align_iterations,
                "pyramid_levels": pyramid_levels,
                "channel_order": ["T1", "FLAIR"],
            },
        }
    )
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Initial alignment checkpoint exists; use resume=True: {manifest_path}")
        old = json.loads(manifest_path.read_text(encoding="utf-8"))
        if old.get("run_signature") != signature:
            raise ValueError("Cannot resume initial atlas alignment: inputs or settings differ")
        cached = _valid_cached(old, signature, root, images)
        if cached is not None:
            return cached
    elif resume:
        raise FileNotFoundError(f"Cannot resume; no initial alignment manifest exists: {manifest_path}")

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage": "plus-lesion-initial-atlas-alignment",
        "manifest_path": str(manifest_path),
        "run_signature": signature,
        "status": "running",
        "spatial_convention": "RAS",
        "source_order": [
            "preprocess.initialCentroids",
            "preprocess.initialMemberships",
            "atlas.initShapeRegistration",
            "atlas.setTransform(preprocess.getTransform())",
            "atlas.updateShapeRegistration",
            "atlas.registerShapesPyramid",
        ],
        "inputs": {
            "cropped_modalities": {"path": str(cropped_path), "sha256": _sha256(cropped_path)},
            "crop_manifest": {"path": str(crop_manifest), "sha256": _sha256(crop_manifest)},
            "atlas_descriptor": {"path": str(descriptor), "sha256": _sha256(descriptor)},
            "topology_atlas": {"path": str(atlas_topology.path), "sha256": _sha256(atlas_topology.path)},
        },
        "settings": {
            "init_align_iterations": init_align_iterations,
            "pyramid_levels": pyramid_levels,
            "channel_order": ["T1", "FLAIR"],
            "inclusive_crop_bounds": list(inclusive_bounds),
            "half_open_membership_bounds": [list(item) for item in half_open_bounds],
            "padding": padding,
            "image_resolution_mm": list(_spacing(images)),
            "template_resolution_mm": list(_spacing(atlas_topology.image)),
        },
        "equivalence_note": "Source-ordered selected initial-alignment composition; full LesionTOADS and real-image workflow validation remain open.",
    }
    _write_json(manifest_path, manifest)

    alignment = initialize_lesion_atlas_alignment_from_descriptor(
        str(descriptor),
        channels_first,
        transform,
        crop_bounds=half_open_bounds,
        pad_border=padding,
        image_resolutions_mm=_spacing(images),
        template_resolutions_mm=_spacing(atlas_topology.image),
        init_align_iterations=init_align_iterations,
        pyramid_levels=pyramid_levels,
    )
    memberships = alignment.alignment.memberships.memberships
    reference = nib.Nifti1Image(
        np.asarray(images.dataobj[..., 0], dtype=np.float32), images.affine, header=images.header
    )
    memberships_image = nifti_ras_on_reference(reference, memberships, dtype=np.float32)
    memberships_path = save_nifti_ras(
        memberships_image, root / "initial_memberships.nii.gz", dtype=np.float32
    )
    summary_path = root / "initial_alignment.json"
    summary: dict[str, Any] = {
        "spatial_convention": "RAS",
        "class_names": [item.name for item in alignment.priors.structures],
        "class_labels": [int(item.label) for item in alignment.priors.structures],
        "modalities": list(alignment.priors.modalities),
        "call_order": list(alignment.alignment.call_order),
        "initial_centroids": alignment.alignment.initial_centroids.tolist(),
        "registration_transform": [float(value) for value in alignment.alignment.registration.transform],
        "preprocessing_transform": [float(value) for value in alignment.alignment.preprocessing_transform],
        "crop_bounds_half_open": [list(item) for item in half_open_bounds],
        "pad_border": padding,
        "membership_shape_xyz_classes": list(memberships.shape),
        "initial_registration_iterations": init_align_iterations,
        "pyramid_levels": pyramid_levels,
        "equivalence_note": manifest["equivalence_note"],
    }
    _write_json(summary_path, summary)
    manifest["status"] = "complete"
    manifest["outputs"] = {
        "memberships": {"path": str(memberships_path), "sha256": _sha256(memberships_path)},
        "summary": {"path": str(summary_path), "sha256": _sha256(summary_path)},
    }
    _write_json(manifest_path, manifest)
    final_transform = tuple(float(value) for value in alignment.alignment.registration.transform)
    return PlusLesionAtlasInitResult(manifest_path, memberships_path, summary_path, final_transform, False)
