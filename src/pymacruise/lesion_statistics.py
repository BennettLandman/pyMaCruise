"""Selected LesionTOADS statistical updates used during membership iteration."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def compute_constrained_lesion_toads_centroids(
    images: ArrayLike,
    memberships: ArrayLike,
    segmentation: ArrayLike,
    centroids: ArrayLike,
    *,
    fields: ArrayLike | None = None,
    relaxation: float = 0.45,
) -> NDArray[np.float32]:
    """Update centroids with selected ``computeConstrainedCentroids`` semantics.

    For each modality and current hard class, the method accumulates a
    membership-squared weighted mean in source x/y/z traversal order. When
    fields are present, it uses ``field * image`` before accumulation. The
    result is constrained by the selected source's 45% previous/next-centroid
    bounds. Missing classes start from a zero sample mean before those same
    bounds are applied. The source assumes an ascending centroid order; this
    helper preserves its arithmetic for any finite input order.
    """
    values = np.asarray(images, dtype=np.float32)
    weights = np.asarray(memberships, dtype=np.float32)
    labels = np.asarray(segmentation)
    initial = np.asarray(centroids, dtype=np.float32)
    field_values = None if fields is None else np.asarray(fields, dtype=np.float32)
    if values.ndim != 4 or any(size < 1 for size in values.shape[:3]) or values.shape[3] < 1:
        raise ValueError("images must have shape (x, y, z, channels)")
    spatial, channels = values.shape[:3], values.shape[3]
    if weights.ndim != 4 or weights.shape[:3] != spatial or weights.shape[3] < 1:
        raise ValueError("memberships must match the image grid and have at least one class")
    classes = weights.shape[3]
    if labels.shape != spatial or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer grid matching the images")
    if initial.shape != (channels, classes):
        raise ValueError("centroids must have shape (channels, classes)")
    if field_values is not None and field_values.shape != values.shape:
        raise ValueError("fields must match the image grid and channel count")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation labels must index the membership classes")
    if not all(np.isfinite(array).all() for array in (values, weights, initial)):
        raise ValueError("images, memberships, and centroids must be finite")
    if field_values is not None and not np.isfinite(field_values).all():
        raise ValueError("fields must be finite")
    factor = np.float32(relaxation)
    if not np.isfinite(factor) or factor < np.float32(0.0):
        raise ValueError("relaxation must be finite and nonnegative")

    result = initial.copy()
    for channel in range(channels):
        previous = np.float32(-np.inf)
        for class_index in range(classes):
            numerator = np.float32(0.0)
            denominator = np.float32(0.0)
            for x in range(spatial[0]):
                for y in range(spatial[1]):
                    for z in range(spatial[2]):
                        if labels[x, y, z] != class_index:
                            continue
                        membership = weights[x, y, z, class_index]
                        square = np.float32(membership * membership)
                        sample = values[x, y, z, channel]
                        if field_values is not None:
                            sample = np.float32(field_values[x, y, z, channel] * sample)
                        numerator = np.float32(numerator + np.float32(square * sample))
                        denominator = np.float32(denominator + square)
            next_centroid = (
                np.float32(numerator / denominator)
                if denominator > np.float32(0.0)
                else np.float32(0.0)
            )
            current = result[channel, class_index]
            lower = np.float32(
                current - np.float32(factor * np.float32(current - previous))
            )
            if next_centroid < lower:
                next_centroid = lower
            elif class_index + 1 < classes:
                upper = np.float32(
                    current
                    + np.float32(
                        factor * np.float32(result[channel, class_index + 1] - current)
                    )
                )
                if next_centroid > upper:
                    next_centroid = upper
            previous = current
            result[channel, class_index] = next_centroid
    return result
