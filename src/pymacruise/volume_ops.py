"""Geometry-aware voxel-wise operations for medical image volumes."""

from __future__ import annotations

import warnings
from collections.abc import Callable
from typing import TypeAlias

import nibabel as nib
import numpy as np
import numpy.typing as npt

VolumeArray: TypeAlias = npt.NDArray[np.generic]
Operation: TypeAlias = Callable[[VolumeArray, VolumeArray], VolumeArray]


class VolumeGeometryError(ValueError):
    """Raised when images cannot be combined without resampling."""


_OPERATIONS: dict[str, Operation] = {
    "ADD": np.add,
    "SUBTRACT": np.subtract,
    "MULTIPLY": np.multiply,
    "DIVIDE": np.divide,
    "MIN": np.minimum,
    "MAX": np.maximum,
}


def _reorient_to_ras(
    image: nib.spatialimages.SpatialImage, name: str
) -> nib.spatialimages.SpatialImage:
    """Return an image with voxel axes ordered Right-Anterior-Superior."""
    if image.affine is None:
        raise VolumeGeometryError(
            f"{name} has no spatial affine; cannot reorient to RAS"
        )

    orientation = nib.aff2axcodes(image.affine)
    if any(axis is None for axis in orientation):
        raise VolumeGeometryError(
            f"{name} has an indeterminate axis orientation; cannot reorient to RAS"
        )
    if orientation == ("R", "A", "S"):
        return image

    warnings.warn(
        f"{name} input orientation is {orientation}; reorienting the image to RAS.",
        UserWarning,
        stacklevel=3,
    )
    return nib.as_closest_canonical(image)


def calculate_volume(
    image_a: nib.spatialimages.SpatialImage,
    image_b: nib.spatialimages.SpatialImage,
    operation: str,
    *,
    affine_atol: float = 1e-5,
) -> nib.spatialimages.SpatialImage:
    """Apply a generic voxel-wise operation to aligned medical volumes.

    Inputs are first axis-reoriented to RAS when needed, with a user-visible
    warning. Arrays are then combined directly: multiplication does not
    normalize or clamp values. Inputs must have identical shapes and physical
    affines after reorientation; registration/resampling is left to an
    explicit step. The output is RAS-oriented and retains image A's image
    class and header while promoting integer inputs to floating point to
    avoid arithmetic overflow.
    """
    op_name = operation.strip().upper()
    try:
        op = _OPERATIONS[op_name]
    except KeyError as exc:
        available = ", ".join(sorted(_OPERATIONS))
        raise ValueError(
            f"Unsupported operation {operation!r}; choose one of: {available}"
        ) from exc

    image_a = _reorient_to_ras(image_a, "image_a")
    image_b = _reorient_to_ras(image_b, "image_b")

    if image_a.shape != image_b.shape:
        raise VolumeGeometryError(
            f"Volume shapes differ: {image_a.shape} versus {image_b.shape}"
        )
    if not np.allclose(image_a.affine, image_b.affine, rtol=0.0, atol=affine_atol):
        raise VolumeGeometryError(
            "Volume affines differ; explicitly resample one image before combining"
        )

    data_a = np.asarray(image_a.dataobj)
    data_b = np.asarray(image_b.dataobj)
    result_dtype = np.result_type(data_a.dtype, data_b.dtype, np.float32)
    data_a = data_a.astype(result_dtype, copy=False)
    data_b = data_b.astype(result_dtype, copy=False)

    with np.errstate(divide="ignore", invalid="ignore"):
        result = op(data_a, data_b)

    header = image_a.header.copy()
    header.set_data_dtype(result.dtype)
    header.set_slope_inter(None, None)
    return image_a.__class__(result, image_a.affine, header)


def calculate_nifti_files(
    input_a: str,
    input_b: str,
    output: str,
    operation: str,
    *,
    affine_atol: float = 1e-5,
) -> None:
    """Read two NIfTI files, apply an operation, and save a RAS-oriented result."""
    if not output.lower().endswith(".nii.gz"):
        raise ValueError("Volume outputs must use the compressed NIfTI suffix .nii.gz")
    image_a = nib.load(input_a)
    image_b = nib.load(input_b)
    result = calculate_volume(image_a, image_b, operation, affine_atol=affine_atol)
    nib.save(result, output)
