"""RAS-normalized compressed NIfTI image I/O."""

from __future__ import annotations

import warnings
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to, resample_to_output


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


def is_axis_aligned_ras_grid(image: nib.spatialimages.SpatialImage) -> bool:
    """Return whether a 3-D image uses positive, axis-aligned RAS voxel axes."""
    if not is_orthogonal_ras_grid(image):
        return False
    spacing = nib.affines.voxel_sizes(image.affine)
    directions = np.asarray(image.affine[:3, :3], dtype=np.float64) / spacing[np.newaxis, :]
    return bool(np.allclose(directions, np.eye(3), rtol=0, atol=1e-6))


def is_orthogonal_ras_grid(image: nib.spatialimages.SpatialImage) -> bool:
    """Return whether a 3-D RAS grid has mutually orthogonal voxel axes.

    Oblique rotation is allowed; sheared axes are not. This distinction lets
    callers preserve a valid oblique acquisition grid without resampling while
    still rejecting geometry that the selected TOADS spacing model cannot
    represent.
    """
    if len(image.shape) != 3 or nib.aff2axcodes(image.affine) != ("R", "A", "S"):
        return False
    spacing = np.asarray(nib.affines.voxel_sizes(image.affine), dtype=np.float64)
    if not np.isfinite(spacing).all() or np.any(spacing <= 0):
        return False
    directions = np.asarray(image.affine[:3, :3], dtype=np.float64) / spacing[np.newaxis, :]
    return bool(np.allclose(directions.T @ directions, np.eye(3), rtol=0, atol=1e-6))


def regrid_to_orthogonal_ras(
    image: nib.spatialimages.SpatialImage,
    *,
    reference: nib.spatialimages.SpatialImage | None = None,
    interpolation_order: int = 1,
) -> nib.spatialimages.SpatialImage:
    """Resample a 3-D image onto an axis-aligned RAS grid.

    When ``reference`` is omitted, an orthogonal RAS grid covering ``image``
    is created at the input voxel spacing. Labels and masks should use order 0;
    continuous intensity images should use order 1. An already suitable image
    is returned unchanged when it defines its own reference.
    """
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)) or len(image.shape) != 3:
        raise ValueError("Regridding requires a 3D scalar NIfTI image")
    if interpolation_order not in (0, 1, 3):
        raise ValueError("Regridding interpolation order must be 0, 1, or 3")
    source = normalize_image_to_ras(image, source_name="regrid input")
    if reference is None and is_axis_aligned_ras_grid(source):
        return source
    if reference is None:
        spacing = nib.affines.voxel_sizes(source.affine)
        target = resample_to_output(
            source,
            voxel_sizes=tuple(float(value) for value in spacing),
            order=interpolation_order,
            mode="constant",
            cval=0.0,
        )
        return normalize_image_to_ras(target, source_name="regridded output")

    target = normalize_image_to_ras(reference, source_name="regrid reference")
    if not is_axis_aligned_ras_grid(target):
        raise ValueError("Regridding reference must use an axis-aligned RAS grid")
    if source.shape == target.shape and np.allclose(source.affine, target.affine, rtol=0, atol=1e-6):
        return source
    sampled = resample_from_to(
        source,
        (target.shape, target.affine),
        order=interpolation_order,
        mode="constant",
        cval=0.0,
    )
    data = np.asarray(sampled.dataobj)
    if interpolation_order == 0:
        data = data.astype(source.get_data_dtype(), copy=False)
    header = target.header.copy()
    header.set_data_dtype(data.dtype)
    header.set_slope_inter(None, None)
    output = nib.Nifti1Image(data, target.affine, header=header)
    qform, qcode = target.get_qform(coded=True)
    _, scode = target.get_sform(coded=True)
    output.set_qform(qform if qcode else None, code=int(qcode), strip_shears=True)
    output.set_sform(target.affine, code=int(scode) if scode else 1)
    return output


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


def nifti_ras_on_reference(
    reference: nib.spatialimages.SpatialImage,
    values: np.ndarray,
    *,
    dtype: np.dtype | str | None = None,
) -> nib.spatialimages.SpatialImage:
    """Create a 3-D/4-D RAS NIfTI image on a 3-D reference grid."""
    if not isinstance(reference, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("reference must be a NIfTI image")
    ras_reference = normalize_image_to_ras(reference, source_name="output reference")
    if len(ras_reference.shape) != 3:
        raise ValueError("reference must be a 3-D image")
    data = np.asarray(values)
    if data.ndim not in (3, 4) or data.shape[:3] != ras_reference.shape:
        raise ValueError("output values must be 3-D or 4-D and match the reference grid")
    if dtype is not None:
        data = data.astype(np.dtype(dtype), copy=False)
    header = ras_reference.header.copy()
    header.set_data_dtype(data.dtype)
    header.set_slope_inter(None, None)
    image_type = nib.Nifti2Image if isinstance(ras_reference, nib.Nifti2Image) else nib.Nifti1Image
    output = image_type(data, ras_reference.affine, header=header)
    qform, qcode = ras_reference.get_qform(coded=True)
    _, scode = ras_reference.get_sform(coded=True)
    if qcode:
        output.set_qform(qform, code=int(qcode), strip_shears=True)
    else:
        output.set_qform(None, code=0)
    output.set_sform(ras_reference.affine, code=int(scode) if scode else 1)
    return output


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
