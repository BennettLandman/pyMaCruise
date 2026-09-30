"""Fixed-input port of the dormant LesionTOADS covariance update."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class LesionToadsCovariances:
    """Per-class covariance matrices and the shared lesion covariance."""

    class_covariance: NDArray[np.float32]
    lesion_covariance: NDArray[np.float32]


def compute_lesion_toads_covariances(
    images: ArrayLike,
    memberships: ArrayLike,
    segmentation: ArrayLike,
    class_centroids: ArrayLike,
    lesion_membership: ArrayLike,
    lesion_centroid: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    *,
    fields: ArrayLike | None = None,
) -> LesionToadsCovariances:
    """Reproduce selected ``computCovariance()`` for fixed maps and statistics.

    This method is dormant in the selected wrapper (its call is commented
    out). Its selected-source quirks are retained: only LESION classes route
    samples with strictly greater lesion membership to the shared lesion
    matrix; ties remain class samples; BLACKHOLE has no special route; and the
    shared lesion matrix is normalized at the position of each LESION class
    in class order.
    """
    values = np.asarray(images, dtype=np.float32)
    weights = np.asarray(memberships, dtype=np.float32)
    labels = np.asarray(segmentation)
    centroids = np.asarray(class_centroids, dtype=np.float32)
    lesion_map = np.asarray(lesion_membership, dtype=np.float32)
    lesion_mean = np.asarray(lesion_centroid, dtype=np.float32)
    field_values = None if fields is None else np.asarray(fields, dtype=np.float32)
    types = tuple(object_types)

    if values.ndim != 4 or any(size < 1 for size in values.shape[:3]) or values.shape[3] < 1:
        raise ValueError("images must have shape (x, y, z, channels)")
    spatial, channels = values.shape[:3], values.shape[3]
    if weights.ndim != 4 or weights.shape[:3] != spatial or weights.shape[3] < 1:
        raise ValueError("memberships must match the image grid and have at least one class")
    classes = weights.shape[3]
    if labels.shape != spatial or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer grid matching the images")
    if centroids.shape != (channels, classes):
        raise ValueError("class_centroids must have shape (channels, classes)")
    if lesion_map.shape != spatial or lesion_mean.shape != (channels,):
        raise ValueError("lesion membership and centroid must match the image grid/channels")
    if len(types) != classes or any(not isinstance(item, str) for item in types):
        raise ValueError("object_types must provide one string per membership class")
    if field_values is not None and field_values.shape != values.shape:
        raise ValueError("fields must match the image grid and channel count")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation labels must index the membership classes")
    finite = (values, weights, centroids, lesion_map, lesion_mean)
    if not all(np.isfinite(array).all() for array in finite):
        raise ValueError("images, memberships, centroids, and lesion maps must be finite")
    if field_values is not None and not np.isfinite(field_values).all():
        raise ValueError("fields must be finite")

    result = np.zeros((classes, channels, channels), dtype=np.float32)
    lesion_result = np.zeros((channels, channels), dtype=np.float32)
    epsilon = np.float32(1.0e-30)

    for class_index in range(classes):
        result[class_index].fill(np.float32(0.0))
        denominator = np.float32(0.0)
        lesion_denominator = np.float32(0.0)
        is_lesion = types[class_index].upper() == "LESION"
        for x in range(spatial[0]):
            for y in range(spatial[1]):
                for z in range(spatial[2]):
                    if labels[x, y, z] != class_index:
                        continue
                    membership = weights[x, y, z, class_index]
                    lesion_value = lesion_map[x, y, z]
                    lesion_sample = is_lesion and membership < lesion_value
                    if lesion_sample:
                        denominator_weight = np.float32(
                            np.float32(lesion_value * lesion_value) + epsilon
                        )
                    else:
                        denominator_weight = np.float32(
                            np.float32(membership * membership) + epsilon
                        )
                    for c1 in range(channels):
                        if field_values is None:
                            sample1 = values[x, y, z, c1]
                        else:
                            sample1 = np.float32(
                                field_values[x, y, z, c1] * values[x, y, z, c1]
                            )
                        mean1 = lesion_mean[c1] if lesion_sample else centroids[c1, class_index]
                        diff1 = np.float32(sample1 - mean1)
                        for c2 in range(channels):
                            if field_values is None:
                                sample2 = values[x, y, z, c2]
                            else:
                                sample2 = np.float32(
                                    field_values[x, y, z, c2] * values[x, y, z, c2]
                                )
                            mean2 = lesion_mean[c2] if lesion_sample else centroids[c2, class_index]
                            diff2 = np.float32(sample2 - mean2)
                            term = np.float32(denominator_weight * diff1)
                            term = np.float32(term * diff2)
                            if lesion_sample:
                                lesion_result[c1, c2] = np.float32(
                                    lesion_result[c1, c2] + term
                                )
                            else:
                                result[class_index, c1, c2] = np.float32(
                                    result[class_index, c1, c2] + term
                                )
                    if lesion_sample:
                        lesion_denominator = np.float32(
                            lesion_denominator + denominator_weight
                        )
                    else:
                        denominator = np.float32(denominator + denominator_weight)

        for c1 in range(channels):
            for c2 in range(channels):
                if denominator != np.float32(0.0):
                    result[class_index, c1, c2] = np.float32(
                        result[class_index, c1, c2] / denominator
                    )
                else:
                    result[class_index, c1, c2] = np.float32(0.0)
                if is_lesion:
                    if lesion_denominator != np.float32(0.0):
                        lesion_result[c1, c2] = np.float32(
                            lesion_result[c1, c2] / lesion_denominator
                        )
                    else:
                        lesion_result[c1, c2] = np.float32(0.0)

    return LesionToadsCovariances(result, lesion_result)
