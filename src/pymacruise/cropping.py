"""RAS-preserving crop and inverse-crop operations from the TOADS adapters.

The crop contract follows Yuankai Huo's supplied `CropBasedOnGM` adapter and
the compiled JIST `CubicVolumeCropper`; the source notices and provenance
inventory are retained in the repository's upstream-license records.
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras


@dataclass(frozen=True)
class CropParameters:
    """Inclusive source bounds and geometry for a JIST-style crop."""

    original_shape: tuple[int, int, int]
    bounds: tuple[tuple[int, int], tuple[int, int], tuple[int, int]]
    background: float
    original_affine: np.ndarray

    @property
    def cropped_shape(self) -> tuple[int, int, int]:
        return tuple(high - low + 1 for low, high in self.bounds)

    @property
    def offset(self) -> tuple[int, int, int]:
        return tuple(low for low, _ in self.bounds)


def find_crop_parameters(
    mask: nib.spatialimages.SpatialImage,
    *,
    threshold: float = 0.0,
    edge_width: int = 5,
    background: float = 0.0,
) -> CropParameters:
    """Find threshold bounds, expanding each axis by ``edge_width`` voxels.

    The Java source initializes extrema from the image dimensions, uses a
    strict ``value > threshold`` comparison, and does not clip expanded bounds
    to the image. An empty mask produced nonsensical Java bounds; this port
    raises instead of creating a malformed array.
    """
    image = _prepare_image(mask, "crop mask")
    if len(image.shape) != 3:
        raise ValueError("Cropping requires a 3D scalar image")
    if edge_width < 0:
        raise ValueError("edge_width must be non-negative")
    if not np.isfinite(threshold) or not np.isfinite(background):
        raise ValueError("threshold and background must be finite")
    data = np.asarray(image.dataobj)
    if not np.isfinite(data).all():
        raise ValueError("Crop mask contains NaN or infinite voxel values")
    indices = np.where(data > threshold)
    if indices[0].size == 0:
        raise ValueError("Cannot crop from an empty threshold mask")
    bounds = tuple(
        (int(axis.min()) - edge_width, int(axis.max()) + edge_width)
        for axis in indices
    )
    return CropParameters(
        original_shape=tuple(int(size) for size in image.shape),
        bounds=bounds,  # type: ignore[arg-type]
        background=float(background),
        original_affine=image.affine.copy(),
    )


def crop_image(
    image: nib.spatialimages.SpatialImage,
    parameters: CropParameters,
) -> nib.Nifti1Image:
    """Apply a shared crop; out-of-image padding uses its configured value."""
    source = _prepare_image(image, "crop target")
    if len(source.shape) != 3:
        raise ValueError("Cropping requires a 3D scalar image")
    if source.shape != parameters.original_shape:
        raise ValueError(
            f"Crop target shape {source.shape} differs from mask shape "
            f"{parameters.original_shape}"
        )
    if not np.allclose(source.affine, parameters.original_affine, rtol=0.0, atol=1e-5):
        raise ValueError("Crop target and mask do not share the same physical grid")
    shape = parameters.cropped_shape
    if any(size <= 0 for size in shape):
        raise ValueError(f"Invalid crop shape {shape}")
    output = np.full(shape, parameters.background, dtype=np.asarray(source.dataobj).dtype)
    source_data = np.asarray(source.dataobj)
    source_slices: list[slice] = []
    target_slices: list[slice] = []
    for (low, high), size in zip(parameters.bounds, source.shape, strict=True):
        start = max(0, low)
        stop = min(size, high + 1)
        source_slices.append(slice(start, stop))
        target_slices.append(slice(start - low, stop - low))
    if all(s.stop > s.start for s in source_slices):
        output[tuple(target_slices)] = source_data[tuple(source_slices)]
    affine = source.affine.copy()
    translation = np.eye(4)
    translation[:3, 3] = parameters.offset
    affine = affine @ translation
    header = source.header.copy()
    header.set_data_dtype(output.dtype)
    header.set_slope_inter(None, None)
    return nib.Nifti1Image(output, affine, header=header)


def crop_from_mask(
    mask: nib.spatialimages.SpatialImage,
    target: nib.spatialimages.SpatialImage | None = None,
    *,
    threshold: float = 0.0,
    edge_width: int = 5,
    background: float = 0.0,
) -> tuple[nib.Nifti1Image, nib.Nifti1Image | None, CropParameters]:
    """Crop a mask and optional companion using one mask-derived transform."""
    mask_ras = _prepare_image(mask, "crop mask")
    parameters = find_crop_parameters(
        mask_ras,
        threshold=threshold,
        edge_width=edge_width,
        background=background,
    )
    cropped_mask = crop_image(mask_ras, parameters)
    cropped_target = crop_image(target, parameters) if target is not None else None
    return cropped_mask, cropped_target, parameters


def uncrop_image(
    cropped: nib.spatialimages.SpatialImage,
    parameters: CropParameters,
    *,
    legacy_exclusive_max: bool = True,
) -> nib.Nifti1Image:
    """Restore a crop to its source grid, padding outside with background.

    ``legacy_exclusive_max`` defaults to the compiled JIST implementation's
    ``index < max`` checks. That historical off-by-one omits the final row,
    column, and slice of the inclusive crop bounds. Set it false for a complete
    mathematical inverse when desired.
    """
    source = _prepare_image(cropped, "uncrop source")
    if len(source.shape) != 3:
        raise ValueError("Uncropping requires a 3D scalar image")
    expected = parameters.cropped_shape
    if source.shape != expected:
        raise ValueError(f"Cropped image shape {source.shape} differs from {expected}")
    output = np.full(
        parameters.original_shape,
        parameters.background,
        dtype=np.asarray(source.dataobj).dtype,
    )
    source_data = np.asarray(source.dataobj)
    # Compute paired source/destination ranges. JIST tests destination indices
    # against the exclusive max bound even though its crop dimensions are
    # inclusive, so the last crop voxel is intentionally skipped by default.
    destination_slices: list[slice] = []
    source_slices: list[slice] = []
    for (low, high), size in zip(
        parameters.bounds, parameters.original_shape, strict=True
    ):
        stop_bound = high if legacy_exclusive_max else high + 1
        dest_start = max(0, low)
        dest_stop = min(size, stop_bound)
        destination_slices.append(slice(dest_start, max(dest_start, dest_stop)))
        source_start = dest_start - low
        source_stop = source_start + max(0, dest_stop - dest_start)
        source_slices.append(slice(source_start, source_stop))
    if all(s.stop > s.start for s in destination_slices):
        output[tuple(destination_slices)] = source_data[tuple(source_slices)]
    header = source.header.copy()
    header.set_data_dtype(output.dtype)
    header.set_slope_inter(None, None)
    return nib.Nifti1Image(output, parameters.original_affine, header=header)


def _prepare_image(
    image: nib.spatialimages.SpatialImage,
    source_name: str,
) -> nib.spatialimages.SpatialImage:
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError(f"{source_name} must be a NIfTI image")
    return normalize_image_to_ras(image, source_name=source_name)
