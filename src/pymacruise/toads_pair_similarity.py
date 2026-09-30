"""Selected VandyCRUISE class-pair similarity prelude.

This is deliberately limited to the ``similarity[k][m]`` coefficient
construction in ``MultipleToadSegmentation.computeExactMemberships``. It
does not implement membership updates, priors, or the wider TOADS iteration.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_math import toads_optimized_centroid_similarity


@dataclass(frozen=True)
class ToadsPairSimilarity:
    """Pairwise source coefficient plus distances computed by the excerpt.

    ``lesion_distance`` is NaN when the Java branch does not assign it. For
    mixed pairs it records the calculated value even though the selected Java
    code then uses ``centroid_distance`` to form the coefficient.
    """

    similarity: NDArray[np.float32]
    centroid_distance: NDArray[np.float32]
    lesion_distance: NDArray[np.float32]


def toads_class_pair_similarity(
    centroids: ArrayLike,
    class_standard_deviations: ArrayLike,
    image_variance_weights: ArrayLike,
    lesion_centroid: ArrayLike,
    lesion_standard_deviation: ArrayLike,
    object_types: ArrayLike,
    prior_scale: float = 0.05,
) -> ToadsPairSimilarity:
    """Dispatch OBJECT/LESION pairs then construct source coefficients.

    ``object_types`` contains source byte codes (OBJECT=1, LESION=4), or the
    strings ``"OBJECT"``/``"LESION"``. For OBJECT/OBJECT and mixed pairs, the
    ordinary class-centroid distance is used. LESION/LESION sets that distance
    to zero. In either mixed direction the source also calculates a distance
    from ``lesion_centroid`` and ``lesion_standard_deviation``; that value is
    intentionally *not* consumed by the following coefficient expression.

    All scalar operations use float32 and preserve Java's left-to-right
    accumulation and divide order. Zero-deviation IEEE NaNs/infinities are
    retained like Java float arithmetic.
    """
    centers = np.asarray(centroids, dtype=np.float32)
    deviations = np.asarray(class_standard_deviations, dtype=np.float32)
    weights = np.asarray(image_variance_weights, dtype=np.float32)
    lesion_center = np.asarray(lesion_centroid, dtype=np.float32)
    lesion_std = np.asarray(lesion_standard_deviation, dtype=np.float32)
    raw_types = np.asarray(object_types)
    if centers.ndim != 2 or centers.shape[0] == 0 or centers.shape[1] == 0:
        raise ValueError("centroids must be a nonempty (modality, class) array")
    modalities, classes = centers.shape
    if deviations.shape != centers.shape:
        raise ValueError("class_standard_deviations must match centroids")
    if weights.shape != (modalities,):
        raise ValueError("image_variance_weights must match modality count")
    if lesion_center.shape != (modalities,) or lesion_std.shape != (modalities,):
        raise ValueError("lesion centroid/deviation must match modality count")
    if raw_types.shape != (classes,):
        raise ValueError("object_types must contain one value per class")
    if not all(np.isfinite(a).all() for a in (centers, deviations, weights, lesion_center, lesion_std)):
        raise ValueError("centroid, deviation, and weight inputs must be finite")

    def is_lesion(value: object) -> bool:
        if isinstance(value, str):
            normalized = value.upper()
            if normalized in {"LESION", "4"}:
                return True
            if normalized in {"OBJECT", "1"}:
                return False
            raise ValueError(f"unsupported object type {value!r}; expected OBJECT or LESION")
        integer = int(value)
        if integer == 4:
            return True
        if integer == 1:
            return False
        raise ValueError(f"unsupported object type {value!r}; expected source code 1 or 4")

    lesion_flags = tuple(is_lesion(v) for v in raw_types)
    # Reuse the existing Java-order coefficient implementation for the
    # ordinary class-centroid distance matrix. It also preserves singular
    # priorScale and zero-deviation IEEE behavior.
    ordinary = toads_optimized_centroid_similarity(
        centers, deviations, weights, prior_scale=prior_scale
    )
    similarity = np.array(ordinary, dtype=np.float32, copy=True)
    centroid_distance = np.zeros((classes, classes), dtype=np.float32)
    lesion_distance = np.full((classes, classes), np.nan, dtype=np.float32)
    scale = np.float32(prior_scale)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for k in range(classes):
            for m in range(classes):
                if k == m:
                    similarity[k, m] = np.float32(0.0)
                    centroid_distance[k, m] = np.float32(0.0)
                    continue
                k_lesion, m_lesion = lesion_flags[k], lesion_flags[m]
                if k_lesion and m_lesion:
                    distance = np.float32(0.0)
                    # Java assigns centroidDist=0 and skips lesionDist here.
                    coefficient = _coefficient(distance, scale)
                else:
                    distance = np.float32(0.0)
                    for channel in range(modalities):
                        delta = np.float32(centers[channel, k] - centers[channel, m])
                        std_sum = np.float32(deviations[channel, k] + deviations[channel, m])
                        term = np.float32(weights[channel] * delta)
                        term = np.float32(term * delta)
                        term = np.float32(term / std_sum)
                        term = np.float32(term / std_sum)
                        distance = np.float32(distance + term)

                    if k_lesion:
                        lesion_distance[k, m] = _distance(
                            lesion_center, centers[:, m], lesion_std, deviations[:, m], weights
                        )
                    elif m_lesion:
                        lesion_distance[k, m] = _distance(
                            lesion_center, centers[:, k], lesion_std, deviations[:, k], weights
                        )
                    # Preserve source quirk: mixed pair coefficient consumes
                    # centroidDist, not the just-computed lesionDist.
                    coefficient = np.float32(ordinary[k, m])
                centroid_distance[k, m] = distance
                similarity[k, m] = coefficient

    return ToadsPairSimilarity(similarity, centroid_distance, lesion_distance)


def _coefficient(distance: np.float32, prior_scale: np.float32) -> np.float32:
    factor = np.float32(np.float32(np.float32(1.0) - prior_scale) / prior_scale)
    scaled = np.float32(factor * distance)
    denominator = np.float32(np.float32(1.0) + scaled)
    return np.float32(np.float32(1.0) / denominator)


def _distance(
    lesion_center: NDArray[np.float32],
    class_center: NDArray[np.float32],
    lesion_std: NDArray[np.float32],
    class_std: NDArray[np.float32],
    weights: NDArray[np.float32],
) -> np.float32:
    distance = np.float32(0.0)
    for channel in range(lesion_center.shape[0]):
        delta = np.float32(lesion_center[channel] - class_center[channel])
        std_sum = np.float32(lesion_std[channel] + class_std[channel])
        term = np.float32(weights[channel] * delta)
        term = np.float32(term * delta)
        term = np.float32(term / std_sum)
        term = np.float32(term / std_sum)
        distance = np.float32(distance + term)
    return distance
