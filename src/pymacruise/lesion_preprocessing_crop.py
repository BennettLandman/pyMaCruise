"""Selected LesionTOADS normalized-image crop and padding setup."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_crop import ToadsCropResult, crop_toads_images, find_toads_crop_bounds


@dataclass(frozen=True)
class LesionPreprocessingCrop:
    """Bounds, zero-padded channels, and transform after source crop."""

    crop: ToadsCropResult
    background_ratio: np.float32
    padding: int


def crop_lesion_toads_preprocessing_inputs(
    normalized_images: ArrayLike,
    topology: ArrayLike,
    transform: ArrayLike,
    *,
    image_resolution_mm: tuple[float, float, float],
    template_resolution_mm: tuple[float, float, float],
    orientation_signs: tuple[float, float, float] = (1.0, 1.0, 1.0),
    padding: int = 10,
) -> LesionPreprocessingCrop:
    """Compose selected `findCroppingBoundaries()` and `cropImages()`.

    The active wrapper normalizes first, then constructs
    `LesionToadPreprocess(..., bgLabel=1, bgRatio=0, border=10, ...)` and
    calls `findCroppingBoundaries()` followed by `cropImages()`. For the
    selected static `robustMaxMin=false` and normalized state, the image test
    includes each voxel whose normalized value is greater than zero. Atlas
    support with label greater than one is transformed into image space and
    unions with that intensity support. Bounds are inclusive. Image channels
    and topology must already be in RAS voxel-axis order; the rotation/resolution
    and zero-padding kernel is reused from the source-validated TOADS adapter.

    The input is fixed, already-normalized data: this operation does not
    reproduce `normalizeImages()` or compute the preprocessing transform.
    """
    images = np.asarray(normalized_images, dtype=np.float32)
    atlas = np.asarray(topology)
    if images.ndim != 4 or any(size == 0 for size in images.shape[1:]):
        raise ValueError("normalized_images must be a nonempty (channels, x, y, z) array")
    if atlas.ndim != 3:
        raise ValueError("topology must be a 3D RAS voxel array")
    if images.shape[1:] == (0, 0, 0):
        raise ValueError("normalized_images spatial dimensions must be nonempty")
    if not np.isfinite(images).all():
        raise ValueError("normalized_images must be finite")
    if isinstance(padding, (bool, np.bool_)) or not isinstance(padding, (int, np.integer)) or padding < 0:
        raise ValueError("padding must be a nonnegative integer")
    trans = np.asarray(transform, dtype=np.float32)
    if trans.shape != (6,) or not np.isfinite(trans).all():
        raise ValueError("transform must contain six finite rigid values")

    # In the active wrapper bgRatio is exactly zero and normalized is true.
    # The second source disjunct has the same >0 threshold in this mode, so
    # [0,1] is a neutral placeholder for the reusable bounds helper's ranges.
    bounds = find_toads_crop_bounds(
        images,
        atlas,
        trans,
        image_resolution=image_resolution_mm,
        template_resolution=template_resolution_mm,
        orientation_signs=orientation_signs,
        robust_ranges=np.tile(np.asarray([[0.0, 1.0]], dtype=np.float32), (images.shape[0], 1)),
        background_ratio=0.0,
        background_label=1,
        normalized=True,
    )
    crop = crop_toads_images(
        images,
        bounds,
        int(padding),
        trans,
        image_resolution=image_resolution_mm,
        orientation_signs=orientation_signs,
    )
    return LesionPreprocessingCrop(crop, np.float32(0.0), int(padding))
