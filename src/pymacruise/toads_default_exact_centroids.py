"""Compose the selected default-atlas exact centroid update.

The 11-class TOADS2008 descriptor contains one MASK and ten optimized tissue
classes (two CSF, five GM, three WM). This adapter routes the MASK through its
selected generic squared-membership centroid branch and the remaining classes
through the already source-matched optimized PowerTable kernel. It covers one
fixed-input post-thinning centroid stage, not the statistics refresh or full
TOADS iteration.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_math import toads_optimized_centroid_update
from .toads_object_centroids import toads_object_centroid_update
from .toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES


_OPTIMIZED_CLASS_TO_CATEGORY = (0, 0, 1, 1, 1, 1, 1, 2, 2, 2)


def toads_default_exact_centroid_update(
    image_values: ArrayLike,
    memberships: ArrayLike,
    previous_centroids: ArrayLike,
    category_exponents: ArrayLike,
    *,
    intensity_field: ArrayLike | None = None,
    smoothness: float = 0.1,
) -> NDArray[np.float32]:
    """Update the default 11-class centroid vector in selected source order.

    Arrays use ``(x,y,z,channel)`` images/fields, ``(x,y,z,11)`` exact
    memberships, ``(channel,11)`` centroids and ``(channel,3)`` exponents
    ordered CSF, GM, WM. Source order is channel, class, x, y, z. For this
    descriptor, class zero is MASK and follows the generic unpowered branch;
    classes 1–2, 3–7 and 8–10 use the descriptor's optimized CSF, GM and WM
    PowerTable factors respectively. The selected wrapper's centroid mode is
    damped with smoothness 0.1. All statistics and fields are fixed inputs.
    """
    images = np.asarray(image_values, dtype=np.float32)
    mems = np.asarray(memberships, dtype=np.float32)
    old = np.asarray(previous_centroids, dtype=np.float32)
    exponents = np.asarray(category_exponents, dtype=np.float32)
    if images.ndim != 4 or images.shape[-1] < 1:
        raise ValueError("image_values must have shape (x,y,z,channel)")
    channels = images.shape[-1]
    if mems.shape != images.shape[:3] + (len(DEFAULT_ATLAS_CLASS_TYPES),):
        raise ValueError("memberships must match the selected 11-class default atlas order")
    if old.shape != (channels, 11):
        raise ValueError("previous_centroids must have shape (channel,11)")
    if exponents.shape != (channels, 3):
        raise ValueError("category_exponents must have shape (channel,CSF/GM/WM)")
    expanded = np.empty((channels, 10), dtype=np.float32)
    for index, category in enumerate(_OPTIMIZED_CLASS_TO_CATEGORY):
        expanded[:, index] = exponents[:, category]

    field = None if intensity_field is None else np.asarray(intensity_field, dtype=np.float32)
    if field is not None and field.shape != images.shape:
        raise ValueError("intensity_field must match image_values")
    updated = old.copy()
    mask_result = toads_object_centroid_update(
        images,
        mems[..., :1],
        old[:, :1],
        intensity_field=field,
        smoothness=smoothness,
    )
    optimized_result = toads_optimized_centroid_update(
        images,
        mems[..., 1:],
        old[:, 1:],
        expanded,
        smoothness,
        intensity_field=field,
    )
    updated[:, :1] = mask_result
    updated[:, 1:] = optimized_result
    return updated
