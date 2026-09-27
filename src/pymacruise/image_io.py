"""RAS-normalized compressed NIfTI image I/O."""

from __future__ import annotations

import warnings
from pathlib import Path

import nibabel as nib
import numpy as np


class ImageGeometryError(ValueError):
    """Raised when a NIfTI image has unusable or ambiguous geometry."""


def normalize_image_to_ras(
    image: nib.spatialimages.SpatialImage,
    *,
    source_name: str = "image",
    affine_atol: float = 1e-4,
) -> nib.spatialimages.SpatialImage:
    """Axis-reorient an image to RAS without interpolation or location changes."""
    affine = image.affine
    if affine is None or not np.isfinite(affine).all():
        raise ImageGeometryError(f"{source_name} has no finite spatial affine")
    if abs(np.linalg.det(affine[:3, :3])) < 1e-10:
        raise ImageGeometryError(f"{source_name} has a singular spatial affine")

    if isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        qform, qcode = image.get_qform(coded=True)
        sform, scode = image.get_sform(coded=True)
        if qcode and scode and not np.allclose(qform, sform, rtol=0, atol=affine_atol):
            raise ImageGeometryError(
                f"{source_name} has conflicting coded qform and sform matrices"
            )

    orientation = nib.aff2axcodes(affine)
    if any(axis is None for axis in orientation):
        raise ImageGeometryError(f"{source_name} has indeterminate axis orientation")
    if orientation == ("R", "A", "S"):
        return image

    warnings.warn(
        f"{source_name} input orientation is {orientation}; reorienting the image to RAS.",
        UserWarning,
        stacklevel=2,
    )
    normalized = nib.as_closest_canonical(image, enforce_diag=False)
    if isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        normalized.set_qform(
            normalized.affine if qcode else None,
            code=int(qcode),
            strip_shears=True,
        )
        normalized.set_sform(
            normalized.affine if scode else None,
            code=int(scode),
        )
    if nib.aff2axcodes(normalized.affine) != ("R", "A", "S"):
        raise ImageGeometryError(f"{source_name} could not be normalized to RAS")
    return normalized


def load_nifti_ras(
    path: str | Path,
    *,
    affine_atol: float = 1e-4,
) -> nib.spatialimages.SpatialImage:
    """Load a compressed NIfTI file, warning and reorienting when needed."""
    source = Path(path)
    if not source.name.lower().endswith(".nii.gz"):
        raise ValueError("Persisted input images must use the compressed .nii.gz format")
    try:
        image = nib.load(str(source))
    except (OSError, nib.filebasedimages.ImageFileError) as exc:
        raise ValueError(f"Cannot read NIfTI image {source}: {exc}") from exc
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError(f"Input is not a NIfTI image: {source}")
    return normalize_image_to_ras(
        image,
        source_name=str(source),
        affine_atol=affine_atol,
    )


def save_nifti_ras(
    image: nib.spatialimages.SpatialImage,
    path: str | Path,
    *,
    dtype: np.dtype | str | None = None,
    affine_atol: float = 1e-4,
) -> Path:
    """Save an image as `.nii.gz` in RAS orientation and return its path."""
    destination = Path(path)
    if not destination.name.lower().endswith(".nii.gz"):
        raise ValueError("Persisted output images must use the compressed .nii.gz format")
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("Only NIfTI images can be persisted by save_nifti_ras")
    normalized = normalize_image_to_ras(
        image,
        source_name="output image",
        affine_atol=affine_atol,
    )
    data = np.asanyarray(normalized.dataobj)
    if dtype is not None:
        data = data.astype(np.dtype(dtype), copy=False)
    header = normalized.header.copy()
    header.set_data_dtype(data.dtype)
    header.set_slope_inter(None, None)
    output = nib.Nifti1Image(data, normalized.affine, header=header)
    qform, qcode = normalized.get_qform(coded=True)
    _, scode = normalized.get_sform(coded=True)
    if qcode:
        output.set_qform(qform, code=int(qcode), strip_shears=True)
    else:
        output.set_qform(None, code=0)
    output.set_sform(normalized.affine, code=int(scode) if scode else 1)
    destination.parent.mkdir(parents=True, exist_ok=True)
    nib.save(output, str(destination))
    return destination
