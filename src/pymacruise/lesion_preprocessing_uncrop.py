"""P39 evidence adapter for the selected LesionTOADS uncrop lifecycle.

The selected ``LesionToadPreprocess.uncropImages()`` array and rigid-translation
operations are already implemented by :func:`pymacruise.toads_crop.uncrop_toads_images`.
This thin adapter makes that source reuse explicit; it does not recreate the
kernel or model the Java object's ``cropped`` flag and crop-center fields.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from .toads_crop import uncrop_toads_images


LESION_PREPROCESS_SOURCE_SHA256 = "72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08"
LESION_UNCROP_METHOD_NORMALIZED_SHA256 = "d2a489171b23e1cd3ff85f2ce203a93e1f1055ca3727493b883ee8b7f5237b18"


def uncrop_lesion_toads_preprocessing_inputs(
    cropped_images: ArrayLike,
    original_shape: tuple[int, int, int],
    bounds: tuple[int, int, int, int, int, int],
    padding: int,
    transform: ArrayLike,
    *,
    image_resolution_mm: tuple[float, float, float],
    orientation_signs: tuple[float, float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """Restore a cropped LesionTOADS image array and undo its crop translation.

    Arrays use channel-first ``[channel, x, y, z]`` with voxel axes already in
    RAS order. Bounds are inclusive. The shared Java-matched kernel zero-fills
    the original grid, copies only the clipped padded region, and reverses the
    rotation-scaled transform translation. This adapter intentionally omits
    output classification/serialization and the Java object's internal
    ``xCi/yCi/zCi`` and ``cropped`` state updates.
    """
    return uncrop_toads_images(
        np.asarray(cropped_images, dtype=np.float32),
        original_shape,
        bounds,
        padding,
        np.asarray(transform, dtype=np.float32),
        image_resolution=image_resolution_mm,
        orientation_signs=orientation_signs,
    )
