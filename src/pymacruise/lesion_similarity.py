"""Centroid-similarity coefficients used by selected LesionTOADS memberships.

This is the opening numerical substep of the selected
``LesionToadSegmentation.computeMemberships`` implementation. It does not
compute voxel memberships or their energies.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class LesionToadsSimilarity:
    """Class and lesion similarities, plus preserved black-hole state."""

    class_similarity: NDArray[np.float32]
    lesion_similarity: NDArray[np.float32]
    black_hole_similarity: NDArray[np.float32]


def _similarity(centroid_distance: np.float32, prior_scale: np.float32) -> np.float32:
    # Keep Java float operation order: 1/(1 + ((1-priorScale)/priorScale)*d).
    scale_ratio = np.float32(np.float32(np.float32(1.0) - prior_scale) / prior_scale)
    scaled_distance = np.float32(scale_ratio * centroid_distance)
    denominator = np.float32(np.float32(1.0) + scaled_distance)
    return np.float32(np.float32(1.0) / denominator)


def compute_lesion_toads_similarity(
    class_centroid: ArrayLike,
    class_stddev: ArrayLike,
    lesion_centroid: ArrayLike,
    lesion_stddev: ArrayLike,
    channel_weight: ArrayLike,
    object_types: Sequence[str],
    structure_names: Sequence[str],
    *,
    prior_scale: float,
    t1_channel: int = -1,
    black_hole_similarity: ArrayLike,
    finite_zero_variance: bool = False,
) -> LesionToadsSimilarity:
    """Reproduce the selected source's class and lesion similarity setup.

    Centroid arrays use ``(channel, class)`` order. With T1 available, a
    class pair involving the named ``Sulcal-CSF`` structure uses only that
    T1 centroid and standard deviation, without channel weighting. All other
    class pairs use the source weighted multi-channel expression. Lesion
    similarity always uses all channels and is zero for lesion, black-hole,
    and optimized-lesion object types. The selected Java method's black-hole
    similarity update is commented out, so the input vector is returned as
    an unchanged copy.

    Zero standard deviations and a zero prior scale retain source IEEE-754
    float32 behavior (including NaN and infinity results) by default. The
    explicit ``finite_zero_variance`` option treats an exactly matching
    centroid with zero deviation as distance zero, and a nonmatching centroid
    as infinite distance, producing finite limiting similarities for a
    stabilized workflow.
    """
    centroids = np.asarray(class_centroid, dtype=np.float32)
    deviations = np.asarray(class_stddev, dtype=np.float32)
    lesion_center = np.asarray(lesion_centroid, dtype=np.float32)
    lesion_deviation = np.asarray(lesion_stddev, dtype=np.float32)
    weights = np.asarray(channel_weight, dtype=np.float32)
    black_similarity = np.asarray(black_hole_similarity, dtype=np.float32)
    types = tuple(object_types)
    names = tuple(structure_names)
    scale = np.float32(prior_scale)
    if not isinstance(finite_zero_variance, (bool, np.bool_)):
        raise TypeError("finite_zero_variance must be boolean")

    if centroids.ndim != 2 or centroids.shape[0] < 1 or centroids.shape[1] < 1:
        raise ValueError("class_centroid must have shape (channels, classes)")
    channels, classes = centroids.shape
    if deviations.shape != centroids.shape:
        raise ValueError("class_stddev must match class_centroid")
    if lesion_center.shape != (channels,) or lesion_deviation.shape != (channels,):
        raise ValueError("lesion centroid and standard deviation must match channels")
    if weights.shape != (channels,):
        raise ValueError("channel_weight must contain one value per channel")
    if len(types) != classes or not all(isinstance(v, str) and v for v in types):
        raise ValueError("one nonempty object type is required per class")
    if len(names) != classes or not all(isinstance(v, str) and v for v in names):
        raise ValueError("one nonempty structure name is required per class")
    if black_similarity.shape != (classes,):
        raise ValueError("black_hole_similarity must contain one value per class")
    for name, values in (
        ("class_centroid", centroids),
        ("class_stddev", deviations),
        ("lesion_centroid", lesion_center),
        ("lesion_stddev", lesion_deviation),
        ("channel_weight", weights),
        ("black_hole_similarity", black_similarity),
    ):
        if not np.isfinite(values).all():
            raise ValueError(f"{name} must contain only finite values")
    if np.any(deviations < 0.0) or np.any(lesion_deviation < 0.0):
        raise ValueError("standard deviations must be nonnegative")
    if not isinstance(t1_channel, (int, np.integer)) or not -1 <= t1_channel < channels:
        raise ValueError("t1_channel must be -1 or a valid channel index")
    if not np.isfinite(scale):
        raise ValueError("prior_scale must be finite")

    class_similarity = np.empty((classes, classes), dtype=np.float32)
    lesion_similarity = np.empty(classes, dtype=np.float32)
    csf_flag = t1_channel != -1
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for k in range(classes):
            for m in range(classes):
                if m == k:
                    class_similarity[k, m] = np.float32(0.0)
                    continue
                distance = np.float32(0.0)
                if csf_flag and (names[m] == "Sulcal-CSF" or names[k] == "Sulcal-CSF"):
                    delta = np.float32(centroids[t1_channel, k] - centroids[t1_channel, m])
                    term = np.float32(delta * delta)
                    denominator = np.float32(
                        deviations[t1_channel, k] + deviations[t1_channel, m]
                    )
                    if finite_zero_variance and denominator == np.float32(0.0):
                        distance = np.float32(0.0) if term == np.float32(0.0) else np.float32(np.inf)
                    else:
                        distance = np.float32(np.float32(term / denominator) / denominator)
                else:
                    for channel in range(channels):
                        delta = np.float32(centroids[channel, k] - centroids[channel, m])
                        # Source is channelWeight * delta * delta / stddev / stddev;
                        # grouping delta squared before weighting changes some float32 ULPs.
                        weighted_delta = np.float32(weights[channel] * delta)
                        weighted = np.float32(weighted_delta * delta)
                        denominator = np.float32(
                            deviations[channel, k] + deviations[channel, m]
                        )
                        if finite_zero_variance and denominator == np.float32(0.0):
                            contribution = (
                                np.float32(0.0)
                                if weighted == np.float32(0.0)
                                else np.float32(np.inf)
                            )
                        else:
                            contribution = np.float32(
                                np.float32(weighted / denominator) / denominator
                            )
                        distance = np.float32(distance + contribution)
                class_similarity[k, m] = _similarity(distance, scale)

        excluded = {"lesion", "blackhole", "optimized_lesion"}
        for k in range(classes):
            distance = np.float32(0.0)
            for channel in range(channels):
                delta = np.float32(lesion_center[channel] - centroids[channel, k])
                weighted_delta = np.float32(weights[channel] * delta)
                weighted = np.float32(weighted_delta * delta)
                denominator = np.float32(lesion_deviation[channel] + deviations[channel, k])
                if finite_zero_variance and denominator == np.float32(0.0):
                    contribution = (
                        np.float32(0.0)
                        if weighted == np.float32(0.0)
                        else np.float32(np.inf)
                    )
                else:
                    contribution = np.float32(np.float32(weighted / denominator) / denominator)
                distance = np.float32(distance + contribution)
            lesion_similarity[k] = (
                _similarity(distance, scale) if types[k] not in excluded else np.float32(0.0)
            )

    return LesionToadsSimilarity(
        class_similarity=class_similarity,
        lesion_similarity=lesion_similarity,
        black_hole_similarity=black_similarity.copy(),
    )
