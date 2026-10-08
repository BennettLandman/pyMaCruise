"""One composed post-thinning iteration for the alternate TOADS approx=4 path.

This adapter joins the selected source-order approximate centroid refresh,
assigned-class variance refresh, and four-slot membership sweep. The short
``MedicAlgorithmToads08`` constructor sets ``useApprox=false``; this API is
therefore for the alternate ``MultipleToadSegmentation`` constructor and is
not a claim about the production wrapper's exact-membership path.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .toads_approximate_centroids import update_default_toads_approximate_centroids
from .toads_approximate_variances import compute_default_toads_approximate_variances
from .toads_default_approximate_pass import (
    ToadsDefaultApproximatePass,
    toads_default_approximate_membership_pass,
)


@dataclass(frozen=True)
class ToadsDefaultApproximateIteration:
    """Outputs from one post-thinning sparse TOADS iteration."""

    centroids: np.ndarray
    centroid_numerators: np.ndarray
    centroid_denominators: np.ndarray
    standard_deviations: np.ndarray
    variance_numerators: np.ndarray
    variance_denominators: np.ndarray
    membership: ToadsDefaultApproximatePass


def run_default_toads_approximate_iteration(
    images: Any,
    memberships: Any,
    best_class_indices: Any,
    segmentation: Any,
    centroids: Any,
    intensity_field: Any,
    optimized_exponents: Any,
    modality_weights: Any,
    exponent_by_class: Any,
    class_variance_weights: Any,
    transformed_shapes: Any,
    direct_relation: Any,
    smoothing_priors: Any,
    atlas_shape_priors: Any,
    *,
    prior_scale: float,
    prior_factor: float,
    prior_coefficient: float,
    smoothing: float,
    centroid_smoothness: float = 0.1,
) -> ToadsDefaultApproximateIteration:
    """Refresh approximate centroids/deviations, then update four memberships.

    Inputs are the state after an approximate-path thinning transition. Image
    and field arrays use ``(x,y,z,channel)`` order; sparse membership and best
    class arrays use ``(x,y,z,slot)`` order with four slots; segmentation uses
    ``(x,y,z)``; and class statistics use ``(channel,class)`` order. The
    selected default atlas has 11 classes and this membership pass currently
    requires one image channel.

    The function does not mutate caller arrays. Its three stages follow the
    alternate constructor's centroid -> class variance -> membership order.
    Each operation has selected-Java method evidence in its dedicated tests;
    this composition adapter itself is validated against direct stage calls.
    """
    centroid_update = update_default_toads_approximate_centroids(
        images,
        memberships,
        best_class_indices,
        centroids,
        intensity_field,
        optimized_exponents,
        centroid_smoothness=centroid_smoothness,
    )
    variance_update = compute_default_toads_approximate_variances(
        images,
        memberships,
        best_class_indices,
        segmentation,
        centroid_update.centroids,
        intensity_field=intensity_field,
    )
    membership_update = toads_default_approximate_membership_pass(
        images,
        memberships,
        best_class_indices,
        centroid_update.centroids,
        variance_update.standard_deviations,
        class_variance_weights,
        modality_weights,
        exponent_by_class,
        transformed_shapes,
        direct_relation,
        smoothing_priors,
        atlas_shape_priors,
        prior_scale=prior_scale,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        smoothing=smoothing,
        intensity_field=intensity_field,
    )
    return ToadsDefaultApproximateIteration(
        centroid_update.centroids,
        centroid_update.numerators,
        centroid_update.denominators,
        variance_update.standard_deviations,
        variance_update.numerators,
        variance_update.denominators,
        membership_update,
    )
