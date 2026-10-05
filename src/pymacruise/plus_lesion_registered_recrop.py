"""Restartable file checkpoint for the selected post-registration recrop."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from .image_io import load_nifti_ras, nifti_ras_on_reference, save_nifti_ras
from .lesion_atlas_priors import (
    load_lesion_atlas_priors,
    load_lesion_atlas_shape_priors,
    load_lesion_atlas_topology,
)
from .lesion_initial_atlas_alignment import (
    LesionInitialAtlasAlignment,
    LesionLoadedAtlasInitialAlignment,
)
from .lesion_initial_memberships import LesionInitialMemberships
from .lesion_shapes import LesionRegistrationPyramidResult
from .lesion_preprocessing_crop import crop_lesion_toads_preprocessing_inputs
from .lesion_registered_recrop import run_registered_lesion_recrop


@dataclass(frozen=True)
class PlusLesionRecropResult:
    """Persisted modalities and atlas maps after registration recropping."""

    manifest_path: Path
    modalities_path: Path
    topology_path: Path
    shapes_path: Path
    summary_path: Path
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
    spacing = np.linalg.norm(np.asarray(image.affine[:3, :3], dtype=np.float64), axis=0)
    if not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError("image geometry must have finite positive RAS voxel spacing")
    return tuple(float(value) for value in spacing)


def _save_crop(
    reference: nib.spatialimages.SpatialImage,
    values: np.ndarray,
    bounds: tuple[int, int, int, int, int, int],
    padding: int,
    path: Path,
    *,
    dtype: np.dtype | str,
) -> Path:
    x0, _, y0, _, z0, _ = bounds
    offset = np.asarray((x0 - padding, y0 - padding, z0 - padding, 1.0), dtype=np.float64)
    affine = np.asarray(reference.affine, dtype=np.float64).copy()
    affine[:3, 3] = (reference.affine @ offset)[:3]
    image = nib.Nifti1Image(np.asarray(values, dtype=dtype), affine, header=reference.header)
    _, qcode = reference.get_qform(coded=True)
    _, scode = reference.get_sform(coded=True)
    image.set_qform(affine if qcode else None, code=int(qcode), strip_shears=True)
    image.set_sform(affine, code=int(scode) if scode else 1)
    return save_nifti_ras(image, path, dtype=dtype)


def _cache_valid(
    manifest: dict[str, Any], signature: str, root: Path
) -> PlusLesionRecropResult | None:
    if manifest.get("status") != "complete" or manifest.get("run_signature") != signature:
        return None
    records = manifest.get("outputs", {})
    paths: dict[str, Path] = {}
    for name in ("modalities", "topology", "shapes", "summary"):
        record = records.get(name)
        if not isinstance(record, dict):
            return None
        path = Path(record.get("path", ""))
        if not path.is_relative_to(root) or not path.is_file() or path.is_symlink():
            return None
        if _sha256(path) != record.get("sha256"):
            return None
        paths[name] = path
    try:
        modalities = nib.load(str(paths["modalities"]))
        topology = nib.load(str(paths["topology"]))
        shapes = nib.load(str(paths["shapes"]))
        summary = json.loads(paths["summary"].read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError):
        return None
    if (
        len(modalities.shape) != 4
        or modalities.shape[-1] != 2
        or topology.shape != modalities.shape[:3]
        or len(shapes.shape) != 4
        or shapes.shape[:3] != modalities.shape[:3]
        or shapes.shape[-1] != 13
        or any(nib.aff2axcodes(image.affine) != ("R", "A", "S") for image in (modalities, topology, shapes))
        or not np.isfinite(np.asarray(modalities.dataobj)).all()
        or not np.isfinite(np.asarray(shapes.dataobj)).all()
    ):
        return None
    if summary.get("spatial_convention") != "RAS":
        return None
    return PlusLesionRecropResult(
        Path(manifest["manifest_path"]), paths["modalities"], paths["topology"],
        paths["shapes"], paths["summary"], True,
    )


def run_plus_lesion_registered_recrop_checkpoint(
    normalized_modalities_path: str | Path,
    crop_manifest_path: str | Path,
    initial_alignment_manifest_path: str | Path,
    atlas_descriptor_path: str | Path,
    output_dir: str | Path,
    *,
    resume: bool = False,
) -> PlusLesionRecropResult:
    """Run the source-ordered uncrop, second crop, and atlas-grid refresh.

    This file boundary consumes P41f normalized full-grid T1/FLAIR and the
    manifests from P41f and P41g. It reconstructs the selected P41g alignment
    from those immutable inputs, then applies the independent preprocessing
    transform for the second crop and the optimized atlas transform for
    transformed topology and shape priors.
    """
    if not isinstance(resume, bool):
        raise TypeError("resume must be boolean")
    normalized_path = Path(normalized_modalities_path).expanduser().resolve()
    crop_manifest_path = Path(crop_manifest_path).expanduser().resolve()
    align_manifest_path = Path(initial_alignment_manifest_path).expanduser().resolve()
    descriptor = Path(atlas_descriptor_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    for role, path in (
        ("normalized modalities", normalized_path),
        ("crop manifest", crop_manifest_path),
        ("initial alignment manifest", align_manifest_path),
        ("atlas descriptor", descriptor),
    ):
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Missing regular {role}: {path}")
    if not normalized_path.name.lower().endswith(".nii.gz"):
        raise ValueError("normalized modalities must be .nii.gz")
    crop_data = json.loads(crop_manifest_path.read_text(encoding="utf-8"))
    align_data = json.loads(align_manifest_path.read_text(encoding="utf-8"))
    if crop_data.get("status") != "complete" or align_data.get("status") != "complete":
        raise ValueError("crop and initial alignment manifests must be complete")
    normalized_record = crop_data.get("outputs", {}).get("normalized_modalities", {})
    if (
        Path(normalized_record.get("path", "")).resolve() != normalized_path
        or _sha256(normalized_path) != normalized_record.get("sha256")
    ):
        raise ValueError("normalized modalities do not match the crop manifest")
    align_inputs = align_data.get("inputs", {})
    if (
        align_inputs.get("crop_manifest", {}).get("path") != str(crop_manifest_path)
        or align_inputs.get("crop_manifest", {}).get("sha256") != _sha256(crop_manifest_path)
    ):
        raise ValueError("initial alignment does not reference this crop manifest")
    cropped_path = Path(crop_data["outputs"]["cropped_modalities"]["path"]).resolve()
    if _sha256(cropped_path) != crop_data["outputs"]["cropped_modalities"]["sha256"]:
        raise ValueError("cropped modalities no longer match the crop manifest")
    if (
        align_inputs.get("cropped_modalities", {}).get("path") != str(cropped_path)
        or align_inputs.get("cropped_modalities", {}).get("sha256") != _sha256(cropped_path)
        or align_inputs.get("atlas_descriptor", {}).get("path") != str(descriptor)
        or align_inputs.get("atlas_descriptor", {}).get("sha256") != _sha256(descriptor)
    ):
        raise ValueError("initial alignment does not match the crop image or atlas descriptor")
    first_topology_record = crop_data.get("inputs", {}).get("topology", {})
    first_topology_path = Path(first_topology_record.get("path", "")).resolve()
    if (
        not first_topology_path.is_file()
        or first_topology_path.is_symlink()
        or _sha256(first_topology_path) != first_topology_record.get("sha256")
    ):
        raise ValueError("initial topology input is missing or no longer matches the crop manifest")
    alignment_summary_record = align_data.get("outputs", {}).get("summary", {})
    alignment_summary_path = Path(alignment_summary_record.get("path", "")).resolve()
    if (
        not alignment_summary_path.is_file()
        or alignment_summary_path.is_symlink()
        or _sha256(alignment_summary_path) != alignment_summary_record.get("sha256")
    ):
        raise ValueError("initial alignment summary is missing or does not match its manifest")
    alignment_summary = json.loads(alignment_summary_path.read_text(encoding="utf-8"))
    membership_record = align_data.get("outputs", {}).get("memberships", {})
    membership_path = Path(membership_record.get("path", "")).resolve()
    if (
        not membership_path.is_file()
        or membership_path.is_symlink()
        or _sha256(membership_path) != membership_record.get("sha256")
    ):
        raise ValueError("initial membership output is missing or no longer matches its manifest")
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "lesion-registered-recrops.json"
    signature = _signature({
        "stage": "plus-lesion-registered-recrops",
        "inputs": {
            "normalized_modalities": {"path": str(normalized_path), "sha256": _sha256(normalized_path)},
            "crop_manifest": {"path": str(crop_manifest_path), "sha256": _sha256(crop_manifest_path)},
            "initial_alignment_manifest": {"path": str(align_manifest_path), "sha256": _sha256(align_manifest_path)},
            "atlas_descriptor": {"path": str(descriptor), "sha256": _sha256(descriptor)},
        },
    })
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Registered recrop checkpoint exists; use resume=True: {manifest_path}")
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if existing.get("run_signature") != signature:
            raise ValueError("Cannot resume registered recrop: inputs differ")
        cached = _cache_valid(existing, signature, root)
        if cached is not None:
            return cached
    elif resume:
        raise FileNotFoundError(f"Cannot resume; no registered recrop manifest exists: {manifest_path}")

    normalized = load_nifti_ras(normalized_path)
    cropped_image = load_nifti_ras(cropped_path)
    memberships_image = load_nifti_ras(membership_path)
    topology = load_lesion_atlas_topology(descriptor)
    first_topology = load_nifti_ras(first_topology_path)
    modalities_last = np.asarray(normalized.dataobj, dtype=np.float32)
    cropped_modalities = np.asarray(cropped_image.dataobj, dtype=np.float32)
    if len(normalized.shape) != 4 or normalized.shape[-1] != 2:
        raise ValueError("normalized modalities must be 4-D channel-last T1 then FLAIR")
    if not np.isfinite(modalities_last).all() or not np.isfinite(cropped_modalities).all():
        raise ValueError("normalized and cropped modalities must be finite")
    memberships_values = np.asarray(memberships_image.dataobj, dtype=np.float32)
    if (
        memberships_values.ndim != 4
        or memberships_values.shape[-1] != 13
        or memberships_values.shape[:3] != cropped_image.shape[:3]
        or nib.aff2axcodes(memberships_image.affine) != ("R", "A", "S")
        or not np.allclose(memberships_image.affine, cropped_image.affine, rtol=0, atol=1e-4)
        or not np.isfinite(memberships_values).all()
    ):
        raise ValueError("initial memberships must be finite 13-class RAS data on the cropped image grid")
    try:
        padding = int(crop_data["settings"]["padding"])
        bounds_values = tuple(int(value) for value in crop_data["crop"]["bounds"])
        preprocessing_transform = np.asarray(crop_data["crop"]["post_crop_transform"], dtype=np.float32)
        iterations = int(align_data["settings"]["init_align_iterations"])
        levels = int(align_data["settings"]["pyramid_levels"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("upstream manifests are missing selected crop/alignment settings") from exc
    if len(bounds_values) != 6 or preprocessing_transform.shape != (6,):
        raise ValueError("upstream crop manifest has invalid geometry")
    original_dimensions = tuple(int(size) for size in normalized.shape[:3])
    topology_values = np.asarray(first_topology.dataobj)
    spacing = _spacing(normalized)
    template_spacing = _spacing(topology.image)
    initial_crop = crop_lesion_toads_preprocessing_inputs(
        modalities_last.transpose(3, 0, 1, 2),
        topology_values,
        np.zeros(6, dtype=np.float32),
        image_resolution_mm=spacing,
        template_resolution_mm=template_spacing,
        padding=padding,
    )
    if initial_crop.crop.bounds != bounds_values or not np.array_equal(
        initial_crop.crop.images.transpose(1, 2, 3, 0).view(np.uint32), cropped_modalities.view(np.uint32)
    ):
        raise ValueError("normalized image and crop manifest do not reproduce the initial crop")
    priors = load_lesion_atlas_priors(descriptor, ("T1", "FLAIR"))
    shape_priors = load_lesion_atlas_shape_priors(descriptor)
    if priors.structures != shape_priors.structures or len(priors.structures) != 13:
        raise ValueError("intensity and shape atlases must share the selected 13-class order")
    registration_transform = np.asarray(alignment_summary["registration_transform"], dtype=np.float32)
    preprocessing_transform = np.asarray(alignment_summary["preprocessing_transform"], dtype=np.float32)
    initial_centroids = np.asarray(alignment_summary["initial_centroids"], dtype=np.float32)
    if (
        registration_transform.shape != (6,)
        or preprocessing_transform.shape != (6,)
        or initial_centroids.shape != (2, 13)
        or not np.isfinite(registration_transform).all()
        or not np.isfinite(preprocessing_transform).all()
        or not np.isfinite(initial_centroids).all()
    ):
        raise ValueError("initial alignment summary has invalid transform or centroid state")
    active_bounds = tuple((padding, size - padding) for size in memberships_values.shape[:3])
    alignment_core = LesionInitialAtlasAlignment(
        initial_centroids=initial_centroids,
        memberships=LesionInitialMemberships(memberships_values, active_bounds),
        registration=LesionRegistrationPyramidResult(registration_transform, ()),
        preprocessing_transform=preprocessing_transform,
        call_order=tuple(alignment_summary["call_order"]),
        initialization_iterations=iterations,
        pyramid_levels=levels,
        registration_scale_factor=2,
        coupling_factor=np.float32(1.0),
    )
    alignment = LesionLoadedAtlasInitialAlignment(priors, shape_priors, alignment_core)

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage": "plus-lesion-registered-recrops",
        "manifest_path": str(manifest_path),
        "run_signature": signature,
        "status": "running",
        "spatial_convention": "RAS",
        "source_order": [
            "preprocess.uncropImages:preprocessingTransform",
            "preprocess.findCroppingBoundaries:preprocessingTransform",
            "preprocess.cropImages:preprocessingTransform",
            "atlas.setImageInfo:recroppedDimensions",
            "atlas.computeTransformedTemplate:optimizedAtlasTransform",
            "preprocess.updateTransformedTemplate",
        ],
        "settings": {"padding": padding, "channel_order": ["T1", "FLAIR"], "classes": 13},
        "equivalence_note": "Source-ordered synthetic post-registration recrop and atlas-map refresh; iterative LesionTOADS and whole-workflow equivalence remain open.",
    }
    _write_json(manifest_path, manifest)
    result = run_registered_lesion_recrop(
        alignment,
        initial_crop,
        topology,
        original_dimensions,
        image_resolutions_mm=spacing,
        template_resolutions_mm=template_spacing,
    )
    recrop = result.preprocessing_crop
    bounds = tuple(int(value) for value in recrop.crop.bounds)
    modalities_output = _save_crop(
        normalized, recrop.crop.images.transpose(1, 2, 3, 0), bounds, padding,
        root / "recropped_modalities_t1_flair.nii.gz", dtype=np.float32,
    )
    modalities_reference = nib.load(str(modalities_output))
    reference = nib.Nifti1Image(
        np.asarray(modalities_reference.dataobj[..., 0], dtype=np.float32),
        modalities_reference.affine,
        header=modalities_reference.header,
    )
    topology_output = save_nifti_ras(
        nifti_ras_on_reference(reference, result.transformed_template, dtype=np.int8),
        root / "transformed_topology.nii.gz", dtype=np.int8,
    )
    shapes_output = save_nifti_ras(
        nifti_ras_on_reference(reference, result.transformed_shapes, dtype=np.float32),
        root / "transformed_shapes_13.nii.gz", dtype=np.float32,
    )
    summary_path = root / "registered_recrop_summary.json"
    summary = {
        "spatial_convention": "RAS",
        "class_names": [item.name for item in alignment.priors.structures],
        "class_labels": [int(item.label) for item in alignment.priors.structures],
        "preprocessing_transform": [float(value) for value in result.preprocessing_crop.crop.transform],
        "atlas_transform": [float(value) for value in result.registration_bounds.transform],
        "recrop_bounds_inclusive": list(bounds),
        "recropped_shape_xyz": list(result.transformed_template.shape),
        "stage_order": list(result.stage_order),
        "equivalence_note": manifest["equivalence_note"],
    }
    _write_json(summary_path, summary)
    manifest["outputs"] = {
        name: {"path": str(path), "sha256": _sha256(path)}
        for name, path in (
            ("modalities", modalities_output), ("topology", topology_output),
            ("shapes", shapes_output), ("summary", summary_path),
        )
    }
    manifest["status"] = "complete"
    _write_json(manifest_path, manifest)
    return PlusLesionRecropResult(
        manifest_path, modalities_output, topology_output, shapes_output, summary_path, False
    )
