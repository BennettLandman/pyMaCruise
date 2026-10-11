"""Selected LesionTOADS statistical updates used during membership iteration."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_math import toads_powertable_lookup


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


def compute_lesion_toads_channel_weights(
    images: ArrayLike,
    memberships: ArrayLike,
    class_centroids: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    outlier_membership: ArrayLike,
    lesion_centroid: ArrayLike,
    black_hole_centroid: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    *,
    outlier_distance: float,
    fields: ArrayLike | None = None,
    norm_type: str = "Diagonal Covariance",
) -> NDArray[np.float32]:
    """Port selected ``computeChannelWeights()`` channel-weight behavior.

    The selected diagonal-covariance path retains the source's squared-error
    OBJECT/OUTLIER terms and its linear residual sum for lesion-related
    classes. MASK and OUTMASK object types do not contribute in this method.
    Other normalization modes return the source's uniform ``1 / nc`` weights.
    """
    values = np.asarray(images, dtype=np.float32)
    weights = np.asarray(memberships, dtype=np.float32)
    centroids = np.asarray(class_centroids, dtype=np.float32)
    lesion_map = np.asarray(lesion_membership, dtype=np.float32)
    black_map = np.asarray(black_hole_membership, dtype=np.float32)
    outlier_map = np.asarray(outlier_membership, dtype=np.float32)
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
    if centroids.shape != (channels, classes):
        raise ValueError("class_centroids must have shape (channels, classes)")
    if any(array.shape != spatial for array in (lesion_map, black_map, outlier_map)):
        raise ValueError("lesion, black-hole, and outlier maps must match the image grid")
    if lesion_mean.shape != (channels,) or black_mean.shape != (channels,):
        raise ValueError("lesion and black-hole centroids must match image channels")
    if field_values is not None and field_values.shape != values.shape:
        raise ValueError("fields must match the image grid and channel count")
    outlier = np.float32(outlier_distance)
    if not np.isfinite(outlier):
        raise ValueError("outlier_distance must be finite")
    finite_arrays = (
        values, weights, centroids, lesion_map, black_map, outlier_map,
        lesion_mean, black_mean,
    )
    if not all(np.isfinite(array).all() for array in finite_arrays):
        raise ValueError("image, membership, and centroid arrays must be finite")
    if field_values is not None and not np.isfinite(field_values).all():
        raise ValueError("fields must be finite")
    if norm_type != "Diagonal Covariance":
        return np.full(channels, np.float32(1.0) / np.float32(channels), dtype=np.float32)

    epsilon = np.float32(1.0e-30)
    infinite = np.float32(1.0e10)
    raw_weights = np.empty(channels, dtype=np.float32)
    weight_sum = np.float32(0.0)
    lesion_types = {"lesion", "optimized_lesion", "blackhole"}
    for channel in range(channels):
        numerator = np.float32(0.0)
        denominator = np.float32(0.0)
        for class_index, kind in enumerate(types):
            for x in range(spatial[0]):
                for y in range(spatial[1]):
                    for z in range(spatial[2]):
                        membership = weights[x, y, z, class_index]
                        sample = values[x, y, z, channel]
                        if field_values is not None:
                            sample = np.float32(field_values[x, y, z, channel] * sample)
                        if kind == "object":
                            difference = np.float32(sample - centroids[channel, class_index])
                            square_weight = np.float32(
                                np.float32(membership * membership) + epsilon
                            )
                            term = np.float32(
                                np.float32(square_weight * difference) * difference
                            )
                            numerator = np.float32(numerator + term)
                            denominator = np.float32(denominator + square_weight)
                        elif kind == "outlier":
                            difference = np.float32(sample - centroids[channel, class_index])
                            outlier_value = outlier_map[x, y, z]
                            complement = np.float32(np.float32(1.0) - outlier_value)
                            inlier_weight = np.float32(
                                np.float32(
                                    np.float32(
                                        np.float32(membership * membership) * complement
                                    )
                                    * complement
                                )
                                + epsilon
                            )
                            inlier_term = np.float32(
                                np.float32(inlier_weight * difference) * difference
                            )
                            outlier_weight = np.float32(
                                np.float32(
                                    np.float32(membership * membership) * outlier_value
                                )
                                * outlier_value
                            )
                            outlier_term = np.float32(
                                np.float32(outlier_weight * outlier) * outlier
                            )
                            numerator = np.float32(
                                numerator + np.float32(inlier_term + outlier_term)
                            )
                            denominator = np.float32(
                                denominator
                                + np.float32(np.float32(membership * membership) + epsilon)
                            )
                        elif kind in lesion_types:
                            difference = np.float32(sample - centroids[channel, class_index])
                            lesion_difference = np.float32(sample - lesion_mean[channel])
                            black_difference = np.float32(sample - black_mean[channel])
                            lesion_value = lesion_map[x, y, z]
                            black_value = black_map[x, y, z]
                            class_weight = np.float32(
                                np.float32(membership * membership) + epsilon
                            )
                            lesion_weight = np.float32(
                                np.float32(lesion_value * lesion_value) + epsilon
                            )
                            black_weight = np.float32(
                                np.float32(black_value * black_value) + epsilon
                            )
                            class_term = np.float32(class_weight * difference)
                            lesion_term = np.float32(lesion_weight * lesion_difference)
                            black_term = np.float32(black_weight * black_difference)
                            numerator_delta = np.float32(
                                np.float32(class_term + lesion_term) + black_term
                            )
                            numerator = np.float32(numerator + numerator_delta)
                            denominator_delta = np.float32(
                                np.float32(class_weight + lesion_weight) + black_weight
                            )
                            denominator = np.float32(denominator + denominator_delta)
        raw = np.float32(denominator / numerator) if numerator > np.float32(0.0) else infinite
        raw_weights[channel] = raw
        weight_sum = np.float32(weight_sum + raw)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.asarray(raw_weights / weight_sum, dtype=np.float32)


def compute_lesion_toads_noise_standard_deviation(
    images: ArrayLike,
    memberships: ArrayLike,
    segmentation: ArrayLike,
    class_centroids: ArrayLike,
    *,
    fields: ArrayLike | None = None,
) -> NDArray[np.float32]:
    """Port selected ``computeNoiseVariance()`` and broadcast per-channel SD.

    The source excludes hard class 0, pools residuals over all remaining
    classes, then assigns each channel's one pooled standard deviation to
    every class. It uses membership-squared weights with an additive ``1e-30``
    floor in both numerator and denominator.
    """
    values = np.asarray(images, dtype=np.float32)
    weights = np.asarray(memberships, dtype=np.float32)
    labels = np.asarray(segmentation)
    centroids = np.asarray(class_centroids, dtype=np.float32)
    field_values = None if fields is None else np.asarray(fields, dtype=np.float32)
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
    if field_values is not None and field_values.shape != values.shape:
        raise ValueError("fields must match the image grid and channel count")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation labels must index the membership classes")
    if not all(np.isfinite(array).all() for array in (values, weights, centroids)):
        raise ValueError("images, memberships, and centroids must be finite")
    if field_values is not None and not np.isfinite(field_values).all():
        raise ValueError("fields must be finite")

    result = np.zeros((channels, classes), dtype=np.float32)
    epsilon = np.float32(1.0e-30)
    for channel in range(channels):
        numerator = np.float32(0.0)
        denominator = np.float32(0.0)
        for class_index in range(1, classes):
            for x in range(spatial[0]):
                for y in range(spatial[1]):
                    for z in range(spatial[2]):
                        if labels[x, y, z] != class_index:
                            continue
                        membership = weights[x, y, z, class_index]
                        sample = values[x, y, z, channel]
                        if field_values is not None:
                            sample = np.float32(field_values[x, y, z, channel] * sample)
                        difference = np.float32(sample - centroids[channel, class_index])
                        square_weight = np.float32(
                            np.float32(membership * membership) + epsilon
                        )
                        term = np.float32(np.float32(square_weight * difference) * difference)
                        numerator = np.float32(numerator + term)
                        denominator = np.float32(denominator + square_weight)
        if denominator > np.float32(0.0):
            quotient = np.float32(numerator / denominator)
            value = np.float32(np.sqrt(np.float64(quotient)))
            result[channel, :] = value
    return result


def compute_lesion_toads_basic_centroids(
    images: ArrayLike,
    memberships: ArrayLike,
    segmentation: ArrayLike,
    centroids: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    *,
    fields: ArrayLike | None = None,
    outlier_membership: ArrayLike | None = None,
    mode: str = "none",
    smoothness: float = 0.0,
    centroid_links: ArrayLike | None = None,
    prior_centroids: ArrayLike | None = None,
    power_exponents: ArrayLike | None = None,
    modalities: tuple[str, ...] | list[str] | None = None,
    zero_threshold: float = 1.0e-10,
) -> NDArray[np.float32]:
    """Update generic LesionTOADS class centroids for supplied hard labels.

    This ports the generic selected ``computeCentroids()`` branch: class
    accumulators use squared memberships and optional field-corrected samples;
    OUTLIER classes first scale class membership by one minus the supplied
    outlier membership map. Optimized CSF/GM/WM classes use the source
    PowerTable residual weight; optimized CSF is gated to the T1_MPRAGE
    modality. ``linked`` combines accumulators for classes sharing a link id,
    ``damped`` blends with the previous centroid, and ``prior`` adds the
    source's ``smoothness * 1e4`` prior weight. The source RAYLEIGH
    second-moment rule is used in ``none``, ``linked``, ``damped`` and fallback
    modes. LESION, BLACKHOLE, and optimized lesion/black-hole categories have
    additional source branches, so this helper rejects them for dedicated
    ports.

    Accumulation follows source channel/class/x/y/z order and explicitly
    rounds intermediate values to float32 to match Java ``float`` arithmetic.
    """
    values = np.asarray(images, dtype=np.float32)
    weights = np.asarray(memberships, dtype=np.float32)
    labels = np.asarray(segmentation)
    initial = np.asarray(centroids, dtype=np.float32)
    field_values = None if fields is None else np.asarray(fields, dtype=np.float32)
    outlier_values = (
        None if outlier_membership is None else np.asarray(outlier_membership, dtype=np.float32)
    )
    types = tuple(item.lower() if isinstance(item, str) else item for item in object_types)
    if values.ndim != 4 or any(size < 1 for size in values.shape[:3]) or values.shape[3] < 1:
        raise ValueError("images must have shape (x, y, z, channels)")
    spatial, channels = values.shape[:3], values.shape[3]
    if weights.ndim != 4 or weights.shape[:3] != spatial or weights.shape[3] < 1:
        raise ValueError("memberships must match the image grid and have at least one class")
    classes = weights.shape[3]
    if len(types) != classes or any(not isinstance(item, str) for item in types):
        raise ValueError("object_types must provide one string per membership class")
    allowed_types = {"object", "mask", "outlier", "outmask", "rayleigh"}
    optimized_types = {"optimized_csf", "optimized_gm", "optimized_wm"}
    allowed_types |= optimized_types
    if any(item not in allowed_types for item in types):
        raise ValueError("object_types include a category requiring a specialized centroid update")
    if "outlier" in types and outlier_values is None:
        raise ValueError("outlier_membership is required when an OUTLIER class is present")
    if labels.shape != spatial or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer grid matching the images")
    if initial.shape != (channels, classes):
        raise ValueError("centroids must have shape (channels, classes)")
    if field_values is not None and field_values.shape != values.shape:
        raise ValueError("fields must match the image grid and channel count")
    if outlier_values is not None and outlier_values.shape != spatial:
        raise ValueError("outlier_membership must match the image grid")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation labels must index the membership classes")
    finite_arrays = (values, weights, initial)
    if not all(np.isfinite(array).all() for array in finite_arrays):
        raise ValueError("images, memberships, and centroids must be finite")
    if field_values is not None and not np.isfinite(field_values).all():
        raise ValueError("fields must be finite")
    if outlier_values is not None and not np.isfinite(outlier_values).all():
        raise ValueError("outlier_membership must be finite")
    blend = np.float32(smoothness)
    if not np.isfinite(blend):
        raise ValueError("smoothness must be finite")

    powers = None if power_exponents is None else np.asarray(power_exponents, dtype=np.float32)
    optimized_present = any(kind in optimized_types for kind in types)
    modality_names = None if modalities is None else tuple(modalities)
    if optimized_present:
        if powers is None or powers.shape != (channels, classes):
            raise ValueError("power_exponents must have shape (channels, classes) for optimized classes")
        if not np.isfinite(powers).all():
            raise ValueError("power_exponents must be finite")
        if modality_names is None or len(modality_names) != channels or any(
            not isinstance(name, str) for name in modality_names
        ):
            raise ValueError("modalities must provide one string per image channel for optimized classes")
        zero = np.float32(zero_threshold)
        if not np.isfinite(zero) or zero < np.float32(0.0):
            raise ValueError("zero_threshold must be finite and nonnegative")
    else:
        zero = np.float32(zero_threshold)

    normalized_mode = mode.lower() if isinstance(mode, str) else ""
    if not normalized_mode:
        raise ValueError("mode must be a nonempty string")
    links = None if centroid_links is None else np.asarray(centroid_links)
    if normalized_mode == "linked":
        if links is None or links.shape != (channels, classes):
            raise ValueError("centroid_links must have shape (channels, classes) in linked mode")
        if not np.issubdtype(links.dtype, np.integer):
            raise ValueError("centroid_links must contain integer link identifiers")
    priors = None if prior_centroids is None else np.asarray(prior_centroids, dtype=np.float32)
    if normalized_mode == "prior":
        if priors is None or priors.shape != (channels, classes):
            raise ValueError("prior_centroids must have shape (channels, classes) in prior mode")
        if not np.isfinite(priors).all():
            raise ValueError("prior_centroids must be finite")

    result = initial.copy()
    rayleigh_types = {index for index, kind in enumerate(types) if kind == "rayleigh"}
    optimized_factors: dict[tuple[int, int], NDArray[np.float32]] = {}
    if optimized_present:
        assert powers is not None and modality_names is not None
        corrected_values = values if field_values is None else np.asarray(field_values * values, dtype=np.float32)
        for channel in range(channels):
            for class_index, kind in enumerate(types):
                uses_table = kind in {"optimized_gm", "optimized_wm"} or (
                    kind == "optimized_csf" and modality_names[channel] == "T1_MPRAGE"
                )
                if not uses_table:
                    continue
                residual = np.asarray(
                    corrected_values[..., channel] - initial[channel, class_index],
                    dtype=np.float32,
                )
                distance = np.asarray(residual * residual, dtype=np.float32)
                factor = np.ones(spatial, dtype=np.float32)
                above_zero = distance > zero
                if np.any(above_zero):
                    exponent = float(powers[channel, class_index])
                    powered = toads_powertable_lookup(
                        distance[above_zero], exponent, cached_exponent=exponent
                    )
                    factor[above_zero] = np.asarray(
                        powered / distance[above_zero].astype(np.float64), dtype=np.float32
                    )
                optimized_factors[(channel, class_index)] = factor
    one = np.float32(1.0)
    two = np.float32(2.0)
    for channel in range(channels):
        numerators = np.zeros(classes, dtype=np.float32)
        denominators = np.zeros(classes, dtype=np.float32)
        for class_index in range(classes):
            for x in range(spatial[0]):
                for y in range(spatial[1]):
                    for z in range(spatial[2]):
                        if labels[x, y, z] != class_index:
                            continue
                        membership = weights[x, y, z, class_index]
                        if types[class_index] == "outlier":
                            assert outlier_values is not None
                            membership = np.float32(
                                membership * np.float32(np.float32(1.0) - outlier_values[x, y, z])
                            )
                        sample = values[x, y, z, channel]
                        if field_values is not None:
                            sample = np.float32(field_values[x, y, z, channel] * sample)
                        square = np.float32(membership * membership)
                        optimized_factor_map = optimized_factors.get((channel, class_index))
                        if optimized_factor_map is not None:
                            factor = optimized_factor_map[x, y, z]
                            weighted = np.float32(
                                np.float32(np.float32(factor * membership) * membership) * sample
                            )
                            weighted_square = np.float32(
                                np.float32(factor * membership) * membership
                            )
                        elif class_index in rayleigh_types:
                            sample_term = np.float32(sample * sample)
                            weighted = np.float32(square * sample_term)
                            weighted_square = square
                        else:
                            weighted = np.float32(square * sample)
                            weighted_square = square
                        numerators[class_index] = np.float32(
                            numerators[class_index] + weighted
                        )
                        denominators[class_index] = np.float32(
                            denominators[class_index] + weighted_square
                        )

        for class_index in range(classes):
            numerator = numerators[class_index]
            denominator = denominators[class_index]
            rayleigh = class_index in rayleigh_types and normalized_mode != "prior"
            if normalized_mode == "linked":
                assert links is not None
                nsum = numerator
                dsum = denominator
                for linked_class in range(classes):
                    if linked_class != class_index and links[channel, class_index] == links[channel, linked_class]:
                        nsum = np.float32(
                            nsum + np.float32(blend * numerators[linked_class])
                        )
                        dsum = np.float32(
                            dsum + np.float32(blend * denominators[linked_class])
                        )
                if dsum > np.float32(0.0):
                    result[channel, class_index] = (
                        np.float32(
                            np.sqrt(np.float64(np.float32(nsum / np.float32(two * dsum))))
                        )
                        if rayleigh else np.float32(nsum / dsum)
                    )
                else:
                    result[channel, class_index] = np.float32(0.0)
            elif normalized_mode == "damped":
                if denominator > np.float32(0.0):
                    if rayleigh:
                        estimate = np.float32(np.sqrt(
                            np.float64(np.float32(numerator / np.float32(two * denominator)))
                        ))
                        retained = np.float32(np.float32(one - blend) * estimate)
                    else:
                        retained = np.float32(
                            np.float32(np.float32(one - blend) * numerator) / denominator
                        )
                    prior = np.float32(blend * result[channel, class_index])
                    result[channel, class_index] = np.float32(retained + prior)
            elif normalized_mode == "prior":
                assert priors is not None
                prior_weight = np.float32(blend * np.float32(1.0e4))
                nsum = np.float32(numerator + np.float32(prior_weight * priors[channel, class_index]))
                dsum = np.float32(denominator + prior_weight)
                result[channel, class_index] = (
                    np.float32(nsum / dsum) if dsum > np.float32(0.0) else np.float32(0.0)
                )
            elif normalized_mode == "none" or normalized_mode not in {"linked", "damped", "prior"}:
                if denominator > np.float32(0.0):
                    if rayleigh:
                        result[channel, class_index] = np.float32(
                            np.sqrt(np.float64(np.float32(numerator / np.float32(two * denominator))))
                        )
                    else:
                        result[channel, class_index] = np.float32(numerator / denominator)
                else:
                    result[channel, class_index] = np.float32(0.0)
    return result


@dataclass(frozen=True)
class LesionToadsAuxiliaryCentroids:
    """Class, lesion, and black-hole centroids from a selected update."""

    class_centroids: NDArray[np.float32]
    lesion_centroids: NDArray[np.float32] | None
    black_hole_centroids: NDArray[np.float32] | None


def compute_lesion_toads_auxiliary_centroids(
    images: ArrayLike,
    memberships: ArrayLike,
    segmentation: ArrayLike,
    centroids: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    lesion_centroids: ArrayLike,
    black_hole_centroids: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    lesion_power_exponents: ArrayLike,
    *,
    fields: ArrayLike | None = None,
    mode: str = "none",
    smoothness: float = 0.0,
    centroid_links: ArrayLike | None = None,
    prior_centroids: ArrayLike | None = None,
    zero_threshold: float = 1.0e-10,
    optimized_power_exponents: ArrayLike | None = None,
) -> LesionToadsAuxiliaryCentroids:
    """Port selected ordinary and optimized auxiliary-centroid routing.

    For hard-labeled LESION/BLACKHOLE classes, the source compares class,
    lesion, and black-hole memberships. Equal maxima route to the class
    membership first, then lesion, then black-hole. Ordinary LESION samples
    receive the source residual-based PowerTable weight. Optimized lesion
    samples instead share a PowerTable factor from winning-membership squared
    times corrected intensity. Class centroids use the selected
    `computeCentroids()` mode; with the source's hard-coded unlinked auxiliary
    centroid branch, lesion estimates are blended 20% toward the sampled mean
    and black-hole estimates are replaced with the sampled mean. The helper
    accepts fixed membership maps and does not estimate them.
    """
    values = np.asarray(images, dtype=np.float32)
    weights = np.asarray(memberships, dtype=np.float32)
    labels = np.asarray(segmentation)
    initial = np.asarray(centroids, dtype=np.float32)
    lesion_map = np.asarray(lesion_membership, dtype=np.float32)
    black_map = np.asarray(black_hole_membership, dtype=np.float32)
    old_lesion = np.asarray(lesion_centroids, dtype=np.float32)
    old_black = np.asarray(black_hole_centroids, dtype=np.float32)
    exponents = np.asarray(lesion_power_exponents, dtype=np.float32)
    optimized_exponents = (
        None if optimized_power_exponents is None
        else np.asarray(optimized_power_exponents, dtype=np.float32)
    )
    field_values = None if fields is None else np.asarray(fields, dtype=np.float32)
    types = tuple(item.lower() if isinstance(item, str) else item for item in object_types)
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
    if len(types) != classes or any(not isinstance(item, str) for item in types):
        raise ValueError("object_types must provide one string per membership class")
    auxiliary_types = {"lesion", "blackhole", "optimized_lesion", "optimized_blackhole"}
    if any(item not in {"object", "mask", "outmask"} | auxiliary_types for item in types):
        raise ValueError("object_types include a category unsupported by auxiliary centroid routing")
    if not (auxiliary_types & set(types)):
        raise ValueError("object_types must include an auxiliary lesion or black-hole category")
    if lesion_map.shape != spatial or black_map.shape != spatial:
        raise ValueError("lesion and black-hole memberships must match the image grid")
    if old_lesion.shape != (channels,) or old_black.shape != (channels,):
        raise ValueError("lesion and black-hole centroids must match the image channels")
    if exponents.shape != (channels,):
        raise ValueError("lesion_power_exponents must have one value per image channel")
    optimized_present = any(
        kind in {"optimized_lesion", "optimized_blackhole"} for kind in types
    )
    if optimized_present:
        if optimized_exponents is None or optimized_exponents.shape != (channels,):
            raise ValueError("optimized_power_exponents must have one value per image channel")
        if not np.isfinite(optimized_exponents).all():
            raise ValueError("optimized_power_exponents must be finite")
    if field_values is not None and field_values.shape != values.shape:
        raise ValueError("fields must match the image grid and channel count")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation labels must index the membership classes")
    finite_arrays = (values, weights, initial, lesion_map, black_map, old_lesion, old_black, exponents)
    if not all(np.isfinite(array).all() for array in finite_arrays):
        raise ValueError("images, maps, centroids, and exponents must be finite")
    if field_values is not None and not np.isfinite(field_values).all():
        raise ValueError("fields must be finite")
    blend = np.float32(smoothness)
    zero = np.float32(zero_threshold)
    if not np.isfinite(blend):
        raise ValueError("smoothness must be finite")
    if not np.isfinite(zero) or zero < np.float32(0.0):
        raise ValueError("zero_threshold must be finite and nonnegative")
    normalized_mode = mode.lower() if isinstance(mode, str) else ""
    if not normalized_mode:
        raise ValueError("mode must be a nonempty string")
    links = None if centroid_links is None else np.asarray(centroid_links)
    if normalized_mode == "linked":
        if links is None or links.shape != (channels, classes):
            raise ValueError("centroid_links must have shape (channels, classes) in linked mode")
        if not np.issubdtype(links.dtype, np.integer):
            raise ValueError("centroid_links must contain integer link identifiers")
    priors = None if prior_centroids is None else np.asarray(prior_centroids, dtype=np.float32)
    if normalized_mode == "prior":
        if priors is None or priors.shape != (channels, classes):
            raise ValueError("prior_centroids must have shape (channels, classes) in prior mode")
        if not np.isfinite(priors).all():
            raise ValueError("prior_centroids must be finite")

    class_result = initial.copy()
    lesion_result = old_lesion.copy()
    black_result = old_black.copy()
    for channel in range(channels):
        numerators = np.zeros(classes, dtype=np.float32)
        denominators = np.zeros(classes, dtype=np.float32)
        lesion_num = np.float32(0.0)
        lesion_den = np.float32(0.0)
        black_num = np.float32(0.0)
        black_den = np.float32(0.0)
        for class_index, kind in enumerate(types):
            for x in range(spatial[0]):
                for y in range(spatial[1]):
                    for z in range(spatial[2]):
                        if labels[x, y, z] != class_index:
                            continue
                        class_membership = weights[x, y, z, class_index]
                        sample = values[x, y, z, channel]
                        if field_values is not None:
                            sample = np.float32(field_values[x, y, z, channel] * sample)
                        if kind in {
                            "lesion", "blackhole", "optimized_lesion", "optimized_blackhole"
                        }:
                            lesion_value = lesion_map[x, y, z]
                            black_value = black_map[x, y, z]
                            winner = max(class_membership, lesion_value, black_value)
                            optimized_kind = kind in {
                                "optimized_lesion", "optimized_blackhole"
                            }
                            factor = np.float32(1.0)
                            if optimized_kind:
                                assert optimized_exponents is not None
                                distance = np.float32(np.float32(winner * winner) * sample)
                                exponent = float(optimized_exponents[channel])
                            elif winner != class_membership and winner == lesion_value:
                                difference = np.float32(sample - old_lesion[channel])
                                distance = np.float32(difference * difference)
                                exponent = float(exponents[channel])
                            else:
                                distance = zero
                                exponent = 0.0
                            if (
                                optimized_kind or (winner != class_membership and winner == lesion_value)
                            ) and distance > zero:
                                powered = toads_powertable_lookup(
                                    distance, exponent, cached_exponent=exponent
                                )
                                factor = np.float32(
                                    np.float64(powered) / np.float64(distance)
                                )
                            if winner == class_membership:
                                weighted = np.float32(
                                    np.float32(np.float32(factor * class_membership)
                                               * class_membership)
                                )
                                numerators[class_index] = np.float32(
                                    numerators[class_index] + np.float32(weighted * sample)
                                )
                                denominators[class_index] = np.float32(
                                    denominators[class_index] + weighted
                                )
                            elif winner == lesion_value:
                                weighted = np.float32(
                                    np.float32(np.float32(factor * lesion_value) * lesion_value)
                                )
                                lesion_num = np.float32(
                                    lesion_num + np.float32(weighted * sample)
                                )
                                lesion_den = np.float32(lesion_den + weighted)
                            else:
                                weighted = np.float32(
                                    np.float32(np.float32(factor * black_value) * black_value)
                                )
                                black_num = np.float32(black_num + np.float32(weighted * sample))
                                black_den = np.float32(black_den + weighted)
                        else:
                            square = np.float32(class_membership * class_membership)
                            term = np.float32(square * sample)
                            numerators[class_index] = np.float32(
                                numerators[class_index] + term
                            )
                            denominators[class_index] = np.float32(
                                denominators[class_index] + square
                            )

        for class_index in range(classes):
            numerator = numerators[class_index]
            denominator = denominators[class_index]
            if normalized_mode == "linked":
                assert links is not None
                nsum = numerator
                dsum = denominator
                for linked_class in range(classes):
                    if linked_class != class_index and links[channel, linked_class] == links[channel, class_index]:
                        nsum = np.float32(nsum + np.float32(blend * numerators[linked_class]))
                        dsum = np.float32(dsum + np.float32(blend * denominators[linked_class]))
                class_result[channel, class_index] = (
                    np.float32(nsum / dsum) if dsum > np.float32(0.0) else np.float32(0.0)
                )
            elif normalized_mode == "damped":
                if denominator > np.float32(0.0):
                    retained = np.float32(
                        np.float32(np.float32(1.0) - blend) * numerator
                    )
                    retained = np.float32(retained / denominator)
                    prior_value = np.float32(blend * class_result[channel, class_index])
                    class_result[channel, class_index] = np.float32(retained + prior_value)
            elif normalized_mode == "prior":
                assert priors is not None
                prior_weight = np.float32(blend * np.float32(1.0e4))
                nsum = np.float32(numerator + np.float32(prior_weight * priors[channel, class_index]))
                dsum = np.float32(denominator + prior_weight)
                class_result[channel, class_index] = (
                    np.float32(nsum / dsum) if dsum > np.float32(0.0) else np.float32(0.0)
                )
            elif denominator > np.float32(0.0):
                class_result[channel, class_index] = np.float32(numerator / denominator)
            elif normalized_mode != "damped":
                class_result[channel, class_index] = np.float32(0.0)

        if lesion_den > np.float32(0.0):
            lesion_mean = np.float32(lesion_num / lesion_den)
            lesion_result[channel] = np.float32(
                np.float32(np.float32(0.2) * lesion_mean)
                + np.float32(np.float32(0.8) * old_lesion[channel])
            )
        else:
            lesion_result[channel] = np.float32(0.0)
        black_result[channel] = (
            np.float32(black_num / black_den)
            if black_den > np.float32(0.0) else np.float32(0.0)
        )
    return LesionToadsAuxiliaryCentroids(class_result, lesion_result, black_result)


def compute_lesion_toads_centroids(
    images: ArrayLike,
    memberships: ArrayLike,
    segmentation: ArrayLike,
    centroids: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    *,
    fields: ArrayLike | None = None,
    outlier_membership: ArrayLike | None = None,
    lesion_membership: ArrayLike | None = None,
    black_hole_membership: ArrayLike | None = None,
    lesion_centroids: ArrayLike | None = None,
    black_hole_centroids: ArrayLike | None = None,
    mode: str = "none",
    smoothness: float = 0.0,
    centroid_links: ArrayLike | None = None,
    prior_centroids: ArrayLike | None = None,
    tissue_power_exponents: ArrayLike | None = None,
    modalities: tuple[str, ...] | list[str] | None = None,
    optimized_csf_all_modalities: bool = False,
    lesion_power_exponents: ArrayLike | None = None,
    optimized_auxiliary_power_exponents: ArrayLike | None = None,
    zero_threshold: float = 1.0e-10,
) -> LesionToadsAuxiliaryCentroids:
    """Update mixed LesionTOADS class and auxiliary centroids in one pass.

    This combines the category-specific branches selected from
    ``LesionToadSegmentation.computeCentroids()`` so linked centroid modes
    operate over one shared class set. It supports generic, OUTLIER,
    RAYLEIGH, optimized CSF/GM/WM, ordinary LESION/BLACKHOLE, and optimized
    LESION/BLACKHOLE classes. Category-specific PowerTable argument sources,
    float32 accumulation order, winner/tie routing, and the source's unlinked
    auxiliary centroid updates are kept explicit. The result's auxiliary
    arrays are ``None`` when no auxiliary class type is present.
    """
    values = np.asarray(images, dtype=np.float32)
    weights = np.asarray(memberships, dtype=np.float32)
    labels = np.asarray(segmentation)
    initial = np.asarray(centroids, dtype=np.float32)
    types = tuple(item.lower() if isinstance(item, str) else item for item in object_types)
    field_values = None if fields is None else np.asarray(fields, dtype=np.float32)
    outlier_values = (
        None if outlier_membership is None
        else np.asarray(outlier_membership, dtype=np.float32)
    )
    lesion_map = None if lesion_membership is None else np.asarray(lesion_membership, dtype=np.float32)
    black_map = None if black_hole_membership is None else np.asarray(black_hole_membership, dtype=np.float32)
    old_lesion = None if lesion_centroids is None else np.asarray(lesion_centroids, dtype=np.float32)
    old_black = None if black_hole_centroids is None else np.asarray(black_hole_centroids, dtype=np.float32)
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
    if len(types) != classes or any(not isinstance(item, str) for item in types):
        raise ValueError("object_types must provide one string per membership class")
    allowed = {
        "object", "mask", "outmask", "outlier", "rayleigh", "optimized_csf",
        "optimized_gm", "optimized_wm", "lesion", "blackhole",
        "optimized_lesion", "optimized_blackhole",
    }
    if any(kind not in allowed for kind in types):
        raise ValueError("object_types contain an unsupported LesionTOADS category")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation labels must index the membership classes")
    if field_values is not None and field_values.shape != values.shape:
        raise ValueError("fields must match the image grid and channel count")
    if "outlier" in types and (outlier_values is None or outlier_values.shape != spatial):
        raise ValueError("OUTLIER classes require an outlier map matching the image grid")
    if outlier_values is not None and outlier_values.shape != spatial:
        raise ValueError("outlier_membership must match the image grid")
    auxiliary_types = {"lesion", "blackhole", "optimized_lesion", "optimized_blackhole"}
    has_auxiliary = any(kind in auxiliary_types for kind in types)
    if has_auxiliary:
        if lesion_map is None or black_map is None or lesion_map.shape != spatial or black_map.shape != spatial:
            raise ValueError("auxiliary class types require lesion and black-hole maps matching the image grid")
        if old_lesion is None or old_black is None or old_lesion.shape != (channels,) or old_black.shape != (channels,):
            raise ValueError("auxiliary class types require one lesion and black-hole centroid per channel")
    if not has_auxiliary and any(item is not None for item in (lesion_map, black_map, old_lesion, old_black)):
        raise ValueError("auxiliary maps and centroids require LESION/BLACKHOLE object types")

    optimized_tissue_types = {"optimized_csf", "optimized_gm", "optimized_wm"}
    has_optimized_tissue = any(kind in optimized_tissue_types for kind in types)
    tissue_powers = None if tissue_power_exponents is None else np.asarray(tissue_power_exponents, dtype=np.float32)
    modality_names = None if modalities is None else tuple(modalities)
    if has_optimized_tissue:
        if tissue_powers is None or tissue_powers.shape != (channels, classes):
            raise ValueError("tissue_power_exponents must have shape (channels, classes)")
        if modality_names is None or len(modality_names) != channels or any(not isinstance(name, str) for name in modality_names):
            raise ValueError("modalities must provide one string per image channel")
    if not isinstance(optimized_csf_all_modalities, bool):
        raise ValueError("optimized_csf_all_modalities must be a boolean")
    ordinary_auxiliary = any(kind in {"lesion", "blackhole"} for kind in types)
    lesion_powers = None if lesion_power_exponents is None else np.asarray(lesion_power_exponents, dtype=np.float32)
    if ordinary_auxiliary and (lesion_powers is None or lesion_powers.shape != (channels,)):
        raise ValueError("lesion_power_exponents must have one value per channel for ordinary auxiliary classes")
    optimized_auxiliary = any(kind in {"optimized_lesion", "optimized_blackhole"} for kind in types)
    optimized_aux_powers = (
        None if optimized_auxiliary_power_exponents is None
        else np.asarray(optimized_auxiliary_power_exponents, dtype=np.float32)
    )
    if optimized_auxiliary and (
        optimized_aux_powers is None or optimized_aux_powers.shape != (channels,)
    ):
        raise ValueError("optimized_auxiliary_power_exponents must have one value per channel")
    finite = [values, weights, initial]
    finite.extend(
        item for item in (
            field_values, outlier_values, lesion_map, black_map, old_lesion,
            old_black, tissue_powers, lesion_powers, optimized_aux_powers,
        ) if item is not None
    )
    if not all(np.isfinite(item).all() for item in finite):
        raise ValueError("images, maps, centroids, and exponents must be finite")
    if not np.isfinite(np.float32(smoothness)):
        raise ValueError("smoothness must be finite")
    blend = np.float32(smoothness)
    zero = np.float32(zero_threshold)
    if not np.isfinite(zero) or zero < np.float32(0.0):
        raise ValueError("zero_threshold must be finite and nonnegative")
    normalized_mode = mode.lower() if isinstance(mode, str) else ""
    if not normalized_mode:
        raise ValueError("mode must be a nonempty string")
    links = None if centroid_links is None else np.asarray(centroid_links)
    if normalized_mode == "linked":
        if links is None or links.shape != (channels, classes) or not np.issubdtype(links.dtype, np.integer):
            raise ValueError("centroid_links must be an integer (channels, classes) map in linked mode")
    priors = None if prior_centroids is None else np.asarray(prior_centroids, dtype=np.float32)
    if normalized_mode == "prior":
        if priors is None or priors.shape != (channels, classes) or not np.isfinite(priors).all():
            raise ValueError("prior_centroids must be finite with shape (channels, classes) in prior mode")

    def _factor_map(distance: NDArray[np.float32], exponent: float) -> NDArray[np.float32]:
        factor = np.ones(spatial, dtype=np.float32)
        selected = distance > zero
        if np.any(selected):
            powered = toads_powertable_lookup(
                distance[selected], exponent, cached_exponent=exponent
            )
            factor[selected] = np.asarray(
                powered / distance[selected].astype(np.float64), dtype=np.float32
            )
        return factor

    result = initial.copy()
    lesion_result = None if old_lesion is None else old_lesion.copy()
    black_result = None if old_black is None else old_black.copy()
    rayleigh_indices = {index for index, kind in enumerate(types) if kind == "rayleigh"}
    for channel in range(channels):
        corrected = values[..., channel] if field_values is None else np.asarray(
            field_values[..., channel] * values[..., channel], dtype=np.float32
        )
        ordinary_lesion_factor = None
        if ordinary_auxiliary:
            residual = np.asarray(corrected - old_lesion[channel], dtype=np.float32)
            lesion_distance = np.asarray(residual * residual, dtype=np.float32)
            ordinary_lesion_factor = _factor_map(lesion_distance, float(lesion_powers[channel]))
        numerators = np.zeros(classes, dtype=np.float32)
        denominators = np.zeros(classes, dtype=np.float32)
        lesion_num = np.float32(0.0)
        lesion_den = np.float32(0.0)
        black_num = np.float32(0.0)
        black_den = np.float32(0.0)
        for class_index, kind in enumerate(types):
            tissue_factor = None
            if kind in optimized_tissue_types:
                # The P41z VandyCRUISE source applies powercsf[c] to every
                # OPTIMIZED_CSF channel. Older selected-source harnesses used
                # a T1_MPRAGE-only policy, retained as the default for those
                # comparisons. The P41z caller opts into the active-source
                # behavior explicitly.
                table_applies = (
                    kind != "optimized_csf"
                    or modality_names[channel] == "T1_MPRAGE"
                    or optimized_csf_all_modalities
                )
                if table_applies:
                    residual = np.asarray(corrected - initial[channel, class_index], dtype=np.float32)
                    distance = np.asarray(residual * residual, dtype=np.float32)
                    tissue_factor = _factor_map(distance, float(tissue_powers[channel, class_index]))
            optimized_aux_factor = None
            if kind in {"optimized_lesion", "optimized_blackhole"}:
                selected_membership = np.maximum(
                    np.maximum(weights[..., class_index], lesion_map), black_map
                )
                distance = np.asarray(
                    np.asarray(selected_membership * selected_membership, dtype=np.float32)
                    * corrected,
                    dtype=np.float32,
                )
                optimized_aux_factor = _factor_map(
                    distance, float(optimized_aux_powers[channel])
                )
            for x in range(spatial[0]):
                for y in range(spatial[1]):
                    for z in range(spatial[2]):
                        if labels[x, y, z] != class_index:
                            continue
                        membership = weights[x, y, z, class_index]
                        sample = corrected[x, y, z]
                        if kind in auxiliary_types:
                            lesion_value = lesion_map[x, y, z]
                            black_value = black_map[x, y, z]
                            winner = max(membership, lesion_value, black_value)
                            factor = np.float32(1.0)
                            optimized_kind = kind in {"optimized_lesion", "optimized_blackhole"}
                            if optimized_kind:
                                factor = optimized_aux_factor[x, y, z]
                            elif winner != membership and winner == lesion_value:
                                factor = ordinary_lesion_factor[x, y, z]
                            if winner == membership:
                                weighted = np.float32(np.float32(factor * membership) * membership)
                                numerators[class_index] = np.float32(
                                    numerators[class_index] + np.float32(weighted * sample)
                                )
                                denominators[class_index] = np.float32(denominators[class_index] + weighted)
                            elif winner == lesion_value:
                                weighted = np.float32(np.float32(factor * lesion_value) * lesion_value)
                                lesion_num = np.float32(lesion_num + np.float32(weighted * sample))
                                lesion_den = np.float32(lesion_den + weighted)
                            else:
                                weighted = np.float32(np.float32(factor * black_value) * black_value)
                                black_num = np.float32(black_num + np.float32(weighted * sample))
                                black_den = np.float32(black_den + weighted)
                            continue
                        if kind == "outlier":
                            membership = np.float32(
                                membership * np.float32(np.float32(1.0) - outlier_values[x, y, z])
                            )
                        square = np.float32(membership * membership)
                        factor = np.float32(1.0) if tissue_factor is None else tissue_factor[x, y, z]
                        weighted = np.float32(np.float32(factor * membership) * membership) if tissue_factor is not None else square
                        sample_term = np.float32(weighted * sample)
                        if class_index in rayleigh_indices:
                            sample_term = np.float32(sample_term * sample)
                        numerators[class_index] = np.float32(
                            numerators[class_index] + sample_term
                        )
                        denominators[class_index] = np.float32(denominators[class_index] + weighted)

        for class_index in range(classes):
            numerator, denominator = numerators[class_index], denominators[class_index]
            if normalized_mode == "linked":
                nsum, dsum = numerator, denominator
                for linked_class in range(classes):
                    if linked_class != class_index and links[channel, class_index] == links[channel, linked_class]:
                        nsum = np.float32(nsum + np.float32(blend * numerators[linked_class]))
                        dsum = np.float32(dsum + np.float32(blend * denominators[linked_class]))
                if dsum > np.float32(0.0):
                    result[channel, class_index] = np.float32(nsum / dsum)
                else:
                    result[channel, class_index] = np.float32(0.0)
            elif normalized_mode == "damped":
                if denominator > np.float32(0.0):
                    retained = np.float32(
                        np.float32(np.float32(np.float32(1.0) - blend) * numerator) / denominator
                    )
                    result[channel, class_index] = np.float32(
                        retained + np.float32(blend * result[channel, class_index])
                    )
            elif normalized_mode == "prior":
                prior_weight = np.float32(blend * np.float32(1.0e4))
                nsum = np.float32(numerator + np.float32(prior_weight * priors[channel, class_index]))
                dsum = np.float32(denominator + prior_weight)
                result[channel, class_index] = np.float32(nsum / dsum) if dsum > np.float32(0.0) else np.float32(0.0)
            elif denominator > np.float32(0.0):
                # The selected Java method applies the Rayleigh square-root
                # conversion only in centroidMode == "none"; linked, damped,
                # prior, and fallback retain their own raw-ratio expressions.
                if class_index in rayleigh_indices and normalized_mode == "none":
                    result[channel, class_index] = np.float32(
                        np.sqrt(np.float64(np.float32(numerator / np.float32(np.float32(2.0) * denominator))))
                    )
                else:
                    result[channel, class_index] = np.float32(numerator / denominator)
            else:
                result[channel, class_index] = np.float32(0.0)

        if has_auxiliary:
            if lesion_den > np.float32(0.0):
                lesion_result[channel] = np.float32(
                    np.float32(np.float32(0.2) * np.float32(lesion_num / lesion_den))
                    + np.float32(np.float32(0.8) * old_lesion[channel])
                )
            else:
                lesion_result[channel] = np.float32(0.0)
            black_result[channel] = (
                np.float32(black_num / black_den)
                if black_den > np.float32(0.0) else np.float32(0.0)
            )
    return LesionToadsAuxiliaryCentroids(result, lesion_result, black_result)
