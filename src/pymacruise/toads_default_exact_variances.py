"""Compose the selected default-atlas exact assigned-class variance update.

All classes in the default TOADS2008 descriptor (MASK and optimized CSF/GM/WM)
share the selected non-LESION branch of ``computeExactVariances``: use full-grid
assigned labels, corrected intensities, the membership-squared plus ``1e-30f``
weight, and source float32 accumulation. The common selected formula is already
implemented by the assigned-class standard-deviation kernel.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES
from .toads_math import toads_optimized_class_standard_deviation


def toads_default_exact_class_variances(
    image_values: ArrayLike,
    memberships: ArrayLike,
    refreshed_centroids: ArrayLike,
    segmentation: ArrayLike,
    *,
    intensity_field: ArrayLike | None = None,
) -> NDArray[np.float32]:
    """Update default 11-class exact standard deviations from fixed state.

    Arrays use ``(x,y,z,channel)`` images/fields, ``(x,y,z,11)`` memberships,
    ``(channel,11)`` refreshed centroids, and an ``(x,y,z)`` label map. The
    selected source scans all voxels for every channel/class, including image
    borders. This composes only ``computeExactVariances`` after centroids have
    been refreshed; it does not estimate centroids or update memberships.
    """
    images = np.asarray(image_values, dtype=np.float32)
    mems = np.asarray(memberships, dtype=np.float32)
    labels = np.asarray(segmentation)
    centers = np.asarray(refreshed_centroids, dtype=np.float32)
    if images.ndim != 4 or images.shape[-1] < 1:
        raise ValueError("image_values must have shape (x,y,z,channel)")
    if mems.shape != images.shape[:3] + (len(DEFAULT_ATLAS_CLASS_TYPES),):
        raise ValueError("memberships must match the selected default 11-class atlas order")
    if labels.shape != images.shape[:3] or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer map matching image_values")
    if centers.shape != (images.shape[-1], 11):
        raise ValueError("refreshed_centroids must have shape (channel,11)")
    return toads_optimized_class_standard_deviation(
        images,
        mems,
        centers,
        labels,
        intensity_field=intensity_field,
    )
