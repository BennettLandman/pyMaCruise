"""Compute centroid similarity in the selected TOADS startup sequence."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_math import toads_optimized_centroid_similarity


@dataclass(frozen=True)
class ToadsInitialCentroidSimilarity:
    """Pairwise class-similarity coefficients from initial statistics."""

    similarity: NDArray[np.float32]
    call_order: tuple[str, ...]


TOADS_INITIAL_SIMILARITY_CALL_ORDER = (
    "segmentation.computeMemberships",
    "segmentation.computeImageVariance",
    "segmentation.computeExactImageVariance",
    "segmentation.computeExactMemberships.centroidSimilarity",
)


def initialize_toads_centroid_similarity(
    *,
    centroids: Any,
    class_standard_deviations: Any,
    image_variance_weights: Any,
    class_types: Any,
    prior_scale: float = 0.05,
) -> ToadsInitialCentroidSimilarity:
    """Run the selected base wrapper's non-lesion similarity prelude.

    VandyCRUISE computes this pairwise coefficient inline in
    ``computeExactMemberships`` immediately after ``computeImageVariance``.
    The base 11-class descriptor uses MASK and optimized tissue classes. This
    adapter preserves all descriptor classes (including MASK) in the matrix;
    lesion-specific pair behavior is outside this base-TOADS slice.
    """
    types = tuple(str(value).lower() for value in class_types)
    if not types or types[0] != "mask" or any(
        value not in {"optimized_csf", "optimized_gm", "optimized_wm"}
        for value in types[1:]
    ):
        raise ValueError("startup similarity expects MASK followed by optimized tissue classes")
    centers = np.asarray(centroids, dtype=np.float32)
    deviations = np.asarray(class_standard_deviations, dtype=np.float32)
    weights = np.asarray(image_variance_weights, dtype=np.float32)
    if centers.ndim != 2 or centers.shape[1] != len(types):
        raise ValueError("centroids must match descriptor class order")
    if deviations.shape != centers.shape:
        raise ValueError("class deviations must match centroid dimensions")
    similarity = toads_optimized_centroid_similarity(
        centers,
        deviations,
        weights,
        prior_scale=prior_scale,
    )
    return ToadsInitialCentroidSimilarity(
        similarity=similarity,
        call_order=TOADS_INITIAL_SIMILARITY_CALL_ORDER,
    )
