"""Selected LesionTOADS statistical updates used during membership iteration."""

from __future__ import annotations

from dataclasses import dataclass

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
        previous = np.float32(-1.0e10)
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


def compute_weighted_lesion_toads_centroids(
    images: ArrayLike,
    memberships: ArrayLike,
    segmentation: ArrayLike,
    centroids: ArrayLike,
    *,
    fields: ArrayLike | None = None,
    factor: float = 0.2,
) -> NDArray[np.float32]:
    """Apply selected ``computeWeightedCentroids(factor)`` semantics.

    The source accumulates membership-squared weighted samples over the
    current hard segmentation, optionally multiplying each image sample by
    its field first, then blends the new mean with the existing centroid.
    If a class has no assigned voxels, the source writes zero without blending.
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
    blend = np.float32(factor)
    if not np.isfinite(blend):
        raise ValueError("factor must be finite")

    result = initial.copy()
    for channel in range(channels):
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
            if denominator > np.float32(0.0):
                retained = np.float32(np.float32(np.float32(1.0) - blend) * numerator)
                estimate = np.float32(retained / denominator)
                prior = np.float32(blend * result[channel, class_index])
                result[channel, class_index] = np.float32(estimate + prior)
            else:
                result[channel, class_index] = np.float32(0.0)
    return result


@dataclass(frozen=True)
class LesionToadsVariances:
    """Class, lesion, and black-hole standard deviations by channel."""

    class_standard_deviation: NDArray[np.float32]
    lesion_standard_deviation: NDArray[np.float32]
    black_hole_standard_deviation: NDArray[np.float32]


def compute_lesion_toads_variances(
    images: ArrayLike,
    memberships: ArrayLike,
    segmentation: ArrayLike,
    class_centroids: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    lesion_centroid: ArrayLike,
    black_hole_centroid: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    *,
    fields: ArrayLike | None = None,
) -> LesionToadsVariances:
    """Port selected ``computeVariances()`` for class/lesion/black-hole maps.

    This deliberately preserves a selected-source detail in the BLACKHOLE
    winner branch: its black-hole residual is weighted by lesion membership,
    not black-hole membership, because that is what the Java method assigns
    to ``lmem`` before updating ``blackHoleStdDev``.
    """
    values = np.asarray(images, dtype=np.float32)
    weights = np.asarray(memberships, dtype=np.float32)
    labels = np.asarray(segmentation)
    centroids = np.asarray(class_centroids, dtype=np.float32)
    lesion_map = np.asarray(lesion_membership, dtype=np.float32)
    black_map = np.asarray(black_hole_membership, dtype=np.float32)
    lesion_mean = np.asarray(lesion_centroid, dtype=np.float32)
    black_mean = np.asarray(black_hole_centroid, dtype=np.float32)
    field_values = None if fields is None else np.asarray(fields, dtype=np.float32)
    types = tuple(object_types)
    if values.ndim != 4 or any(size < 1 for size in values.shape[:3]) or values.shape[3] < 1:
        raise ValueError("images must have shape (x, y, z, channels)")
    spatial, channels = values.shape[:3], values.shape[3]
    if weights.ndim != 4 or weights.shape[:3] != spatial or weights.shape[3] < 1:
        raise ValueError("memberships must match the image grid and have at least one class")
    classes = weights.shape[3]
    if len(types) != classes or any(not isinstance(item, str) for item in types):
        raise ValueError("object_types must provide one string per membership class")
    if labels.shape != spatial or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer grid matching the images")
    if centroids.shape != (channels, classes):
        raise ValueError("class_centroids must have shape (channels, classes)")
    if lesion_map.shape != spatial or black_map.shape != spatial:
        raise ValueError("lesion and black-hole memberships must match the image grid")
    if lesion_mean.shape != (channels,) or black_mean.shape != (channels,):
        raise ValueError("lesion and black-hole centroids must match the image channels")
    if field_values is not None and field_values.shape != values.shape:
        raise ValueError("fields must match the image grid and channel count")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation labels must index the membership classes")
    finite_arrays = (values, weights, centroids, lesion_map, black_map, lesion_mean, black_mean)
    if not all(np.isfinite(array).all() for array in finite_arrays):
        raise ValueError("images, memberships, and centroids must be finite")
    if field_values is not None and not np.isfinite(field_values).all():
        raise ValueError("fields must be finite")

    class_std = np.zeros((channels, classes), dtype=np.float32)
    lesion_std = np.zeros(channels, dtype=np.float32)
    black_std = np.zeros(channels, dtype=np.float32)
    epsilon = np.float32(1.0e-30)
    lesion_types = {"lesion", "optimized_lesion"}
    for channel in range(channels):
        lesion_num = np.float32(0.0)
        lesion_den = np.float32(0.0)
        black_num = np.float32(0.0)
        black_den = np.float32(0.0)
        for class_index in range(classes):
            class_num = np.float32(0.0)
            class_den = np.float32(0.0)
            kind = types[class_index]
            for x in range(spatial[0]):
                for y in range(spatial[1]):
                    for z in range(spatial[2]):
                        if labels[x, y, z] != class_index:
                            continue
                        class_membership = weights[x, y, z, class_index]
                        lesion_value = lesion_map[x, y, z]
                        black_value = black_map[x, y, z]
                        auxiliary = (
                            (kind in lesion_types and class_membership < lesion_value)
                            or (kind == "blackhole" and class_membership < black_value)
                        )
                        sample = values[x, y, z, channel]
                        if field_values is not None:
                            sample = np.float32(field_values[x, y, z, channel] * sample)
                        if auxiliary:
                            if lesion_value >= black_value:
                                selected = lesion_value
                                difference = np.float32(sample - lesion_mean[channel])
                                square_weight = np.float32(
                                    np.float32(selected * selected) + epsilon
                                )
                                term = np.float32(np.float32(square_weight * difference) * difference)
                                lesion_num = np.float32(lesion_num + term)
                                lesion_den = np.float32(lesion_den + square_weight)
                            else:
                                # The source uses lesionMembership here even when black-hole wins.
                                selected = lesion_value
                                difference = np.float32(sample - black_mean[channel])
                                square_weight = np.float32(
                                    np.float32(selected * selected) + epsilon
                                )
                                term = np.float32(np.float32(square_weight * difference) * difference)
                                black_num = np.float32(black_num + term)
                                black_den = np.float32(black_den + square_weight)
                        else:
                            difference = np.float32(sample - centroids[channel, class_index])
                            square_weight = np.float32(
                                np.float32(class_membership * class_membership) + epsilon
                            )
                            term = np.float32(np.float32(square_weight * difference) * difference)
                            class_num = np.float32(class_num + term)
                            class_den = np.float32(class_den + square_weight)
            if class_den > np.float32(0.0):
                quotient = np.float32(class_num / class_den)
                class_std[channel, class_index] = np.float32(np.sqrt(np.float64(quotient)))
            else:
                class_std[channel, class_index] = np.float32(0.0)
        if lesion_den > np.float32(0.0):
            lesion_std[channel] = np.float32(np.sqrt(np.float64(np.float32(lesion_num / lesion_den))))
        else:
            lesion_std[channel] = np.float32(0.0)
        if black_den > np.float32(0.0):
            black_std[channel] = np.float32(np.sqrt(np.float64(np.float32(black_num / black_den))))
        else:
            black_std[channel] = np.float32(0.0)
    return LesionToadsVariances(class_std, lesion_std, black_std)
