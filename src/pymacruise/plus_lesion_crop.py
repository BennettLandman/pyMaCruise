"""Source-ordered LesionTOADS intensity normalization and crop checkpoint."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras, save_nifti_ras
from pymacruise.lesion_preprocessing_crop import crop_lesion_toads_preprocessing_inputs
from pymacruise.lesion_toads_preprocessing import normalize_lesion_toads_images


@dataclass(frozen=True)
class PlusLesionCropResult:
    """Normalized whole-grid and source-padded LesionTOADS modalities."""

    manifest_path: Path
    normalized_path: Path
    cropped_path: Path
    crop_bounds: tuple[int, int, int, int, int, int]
    post_crop_transform: tuple[float, float, float, float, float, float]
    resumed: bool


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _write_manifest(path: Path, manifest: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _same_grid(left: nib.spatialimages.SpatialImage, right: nib.spatialimages.SpatialImage) -> bool:
    return left.shape[:3] == right.shape[:3] and np.allclose(
        left.affine, right.affine, rtol=0, atol=1e-4
    )


def _spacings(image: nib.spatialimages.SpatialImage) -> tuple[float, float, float]:
    spacing = np.linalg.norm(np.asarray(image.affine[:3, :3], dtype=np.float64), axis=0)
    if spacing.shape != (3,) or not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError("image geometry must have finite positive RAS voxel spacing")
    return tuple(float(item) for item in spacing)


def _save_cropped_modalities(
    reference: nib.spatialimages.SpatialImage,
    channel_first_values: np.ndarray,
    bounds: tuple[int, int, int, int, int, int],
    padding: int,
    destination: Path,
) -> Path:
    values = np.asarray(channel_first_values, dtype=np.float32).transpose(1, 2, 3, 0)
    x0, _, y0, _, z0, _ = bounds
    offset = np.asarray((x0 - padding, y0 - padding, z0 - padding, 1.0), dtype=np.float64)
    affine = np.asarray(reference.affine, dtype=np.float64) @ np.diag((1.0, 1.0, 1.0, 1.0))
    affine[:3, 3] = (np.asarray(reference.affine, dtype=np.float64) @ offset)[:3]
    header = reference.header.copy()
    header.set_data_dtype(np.float32)
    header.set_slope_inter(None, None)
    output = nib.Nifti1Image(values, affine, header=header)
    _, qcode = reference.get_qform(coded=True)
    _, scode = reference.get_sform(coded=True)
    output.set_qform(affine if qcode else None, code=int(qcode), strip_shears=True)
    output.set_sform(affine, code=int(scode) if scode else 1)
    return save_nifti_ras(output, destination, dtype=np.float32)


def _valid_result(
    manifest: dict[str, Any], signature: str, root: Path, reference: nib.spatialimages.SpatialImage
) -> PlusLesionCropResult | None:
    if manifest.get("status") != "complete" or manifest.get("run_signature") != signature:
        return None
    records = manifest.get("outputs", {})
    paths: dict[str, Path] = {}
    for name in ("normalized_modalities", "cropped_modalities"):
        record = records.get(name)
        if not isinstance(record, dict):
            return None
        path = Path(record.get("path", ""))
        if not path.is_relative_to(root) or not path.is_file() or path.is_symlink():
            return None
        if _sha256(path) != record.get("sha256"):
            return None
        image = nib.load(str(path))
        if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)) or len(image.shape) != 4:
            return None
        if image.shape[-1] != 2 or nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            return None
        if name == "normalized_modalities" and not _same_grid(image, reference):
            return None
        if not np.isfinite(np.asarray(image.dataobj)).all():
            return None
        paths[name] = path
    try:
        bounds = tuple(int(value) for value in manifest["crop"]["bounds"])
        transform = tuple(float(value) for value in manifest["crop"]["post_crop_transform"])
        padding = int(manifest["settings"]["padding"])
    except (KeyError, TypeError, ValueError):
        return None
    expected = (
        bounds[1] - bounds[0] + 1 + 2 * padding,
        bounds[3] - bounds[2] + 1 + 2 * padding,
        bounds[5] - bounds[4] + 1 + 2 * padding,
        2,
    )
    cropped = nib.load(str(paths["cropped_modalities"]))
    if cropped.shape != expected:
        return None
    return PlusLesionCropResult(
        Path(manifest["manifest_path"]),
        paths["normalized_modalities"],
        paths["cropped_modalities"],
        bounds,
        transform,
        True,
    )


def run_plus_lesion_crop_checkpoint(
    modalities_path: str | Path,
    topology_path: str | Path,
    output_dir: str | Path,
    *,
    padding: int = 10,
    resume: bool = False,
) -> PlusLesionCropResult:
    """Normalize T1/FLAIR, find source crop bounds, and persist the padded crop.

    ``modalities_path`` is the 4-D channel-last output of
    :func:`run_plus_lesion_preprocessing_checkpoint`; its last axis must be
    ``(T1, FLAIR)``. The supplied topology atlas is loaded/reoriented to RAS.
    Since the loaded voxel axes are RAS, the selected axial initial transform
    and orientation signs are zero and ``(+1,+1,+1)`` respectively. The crop
    unions all normalized nonzero image support with transformed atlas labels
    above the source background label, then uses the selected zero padding.
    """
    if isinstance(padding, bool) or not isinstance(padding, int) or padding < 0:
        raise ValueError("padding must be a nonnegative integer")
    if not isinstance(resume, bool):
        raise TypeError("resume must be boolean")
    modalities_file = Path(modalities_path).expanduser().resolve()
    topology_file = Path(topology_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    for role, path in (("modalities", modalities_file), ("topology", topology_file)):
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Missing regular {role} input: {path}")
        if not path.name.lower().endswith(".nii.gz"):
            raise ValueError(f"{role} input must be .nii.gz")

    modalities = load_nifti_ras(modalities_file)
    topology = load_nifti_ras(topology_file)
    if len(modalities.shape) != 4 or modalities.shape[-1] != 2:
        raise ValueError("modalities must be 4-D with channel-last order (T1, FLAIR)")
    if len(topology.shape) != 3:
        raise ValueError("topology atlas must be a 3-D label image")
    values = np.asarray(modalities.dataobj, dtype=np.float32)
    atlas = np.asarray(topology.dataobj)
    if not np.isfinite(values).all() or not np.isfinite(atlas).all():
        raise ValueError("modalities and topology atlas must contain finite values")
    if not np.issubdtype(atlas.dtype, np.number):
        raise ValueError("topology atlas must contain numeric labels")
    channels_first = values.transpose(3, 0, 1, 2)
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "lesion-crop.json"
    signature = _json_hash(
        {
            "stage": "plus-lesion-normalization-crop",
            "inputs": {
                "modalities": {"path": str(modalities_file), "sha256": _sha256(modalities_file)},
                "topology": {"path": str(topology_file), "sha256": _sha256(topology_file)},
            },
            "padding": padding,
            "initial_transform": [0.0] * 6,
            "source_channel_order": ["T1", "FLAIR"],
            "spatial_convention": "RAS",
        }
    )
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Lesion crop checkpoint exists; use resume=True: {manifest_path}")
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if existing.get("run_signature") != signature:
            raise ValueError("Cannot resume lesion crop: inputs or settings differ")
        cached = _valid_result(existing, signature, root, modalities)
        if cached is not None:
            return cached
    elif resume:
        raise FileNotFoundError(f"Cannot resume; no lesion crop manifest exists: {manifest_path}")

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage": "plus-lesion-normalization-crop",
        "manifest_path": str(manifest_path),
        "run_signature": signature,
        "status": "running",
        "spatial_convention": "RAS",
        "source_order": [
            "normalizeImages",
            "findCroppingBoundaries",
            "cropImages",
            "update crop transform translation",
        ],
        "inputs": {
            "modalities": {"path": str(modalities_file), "sha256": _sha256(modalities_file)},
            "topology": {"path": str(topology_file), "sha256": _sha256(topology_file)},
        },
        "settings": {
            "padding": padding,
            "initial_transform": [0.0] * 6,
            "background_label": 1,
            "background_ratio": 0.0,
            "image_resolution_mm": list(_spacings(modalities)),
            "template_resolution_mm": list(_spacings(topology)),
            "orientation_signs": [1.0, 1.0, 1.0],
        },
        "channel_order": ["T1", "FLAIR"],
        "equivalence_note": "Source-ordered synthetic preparation composition; no complete LesionTOADS execution or whole-workflow equivalence claim.",
    }
    _write_manifest(manifest_path, manifest)

    normalized = normalize_lesion_toads_images(channels_first)
    initial_transform = np.zeros(6, dtype=np.float32)
    crop = crop_lesion_toads_preprocessing_inputs(
        normalized.normalized_images,
        atlas,
        initial_transform,
        image_resolution_mm=_spacings(modalities),
        template_resolution_mm=_spacings(topology),
        orientation_signs=(1.0, 1.0, 1.0),
        padding=padding,
    )
    normalized_last = normalized.normalized_images.transpose(1, 2, 3, 0)
    normalized_output = save_nifti_ras(
        nib.Nifti1Image(normalized_last, modalities.affine, header=modalities.header),
        root / "normalized_modalities_t1_flair.nii.gz",
        dtype=np.float32,
    )
    cropped_output = _save_cropped_modalities(
        modalities,
        crop.crop.images,
        crop.crop.bounds,
        padding,
        root / "cropped_modalities_t1_flair.nii.gz",
    )
    manifest["status"] = "complete"
    manifest["crop"] = {
        "bounds": list(crop.crop.bounds),
        "post_crop_transform": [float(value) for value in crop.crop.transform],
        "padded_shape_xyz": list(crop.crop.images.shape[1:]),
    }
    manifest["intensity_bounds"] = {
        "minimum": [float(value) for value in normalized.intensity_minima],
        "maximum": [float(value) for value in normalized.intensity_maxima],
    }
    manifest["outputs"] = {
        "normalized_modalities": {"path": str(normalized_output), "sha256": _sha256(normalized_output)},
        "cropped_modalities": {"path": str(cropped_output), "sha256": _sha256(cropped_output)},
    }
    _write_manifest(manifest_path, manifest)
    return PlusLesionCropResult(
        manifest_path,
        normalized_output,
        cropped_output,
        crop.crop.bounds,
        tuple(float(value) for value in crop.crop.transform),
        False,
    )
