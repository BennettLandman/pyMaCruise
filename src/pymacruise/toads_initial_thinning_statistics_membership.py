"""Continue selected 11-class TOADS startup through statistics and membership.

The caller supplies the already-computed P18cyn first-thinning composition.
This bounded stage refreshes exact centroids and assigned-class deviations
from that thinned segmentation, recomputes centroid similarity, and executes
one subsequent exact-membership pass. Registration and estimation of the
other fixed statistics remain outside this seam.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_default_atlas_membership import (
    ToadsDefaultAtlasMembershipPass,
    toads_default_atlas_membership_pass,
)
from .toads_default_exact_centroids import toads_default_exact_centroid_update
from .toads_default_exact_variances import toads_default_exact_class_variances
from .toads_initial_thinning_composition import ToadsInitialThinningComposition
from .toads_math import toads_optimized_centroid_similarity


@dataclass(frozen=True)
class ToadsInitialStatisticsMembershipContinuation:
    """Selected startup, first thinning, refreshed statistics, and next pass."""

    startup_and_thinning: ToadsInitialThinningComposition
    centroids: NDArray[np.float32]
    class_standard_deviations: NDArray[np.float32]
    centroid_similarity: NDArray[np.float32]
    membership_pass: ToadsDefaultAtlasMembershipPass
    stage_order: tuple[str, ...]


def continue_initial_thinning_through_memberships(
    startup_and_thinning: ToadsInitialThinningComposition,
    *,
    image_values: Any,
    initial_centroids: Any,
    category_exponents: Any,
    transformed_shapes: Any,
    image_variance_weights: Any,
    smoothing_priors: Any,
    atlas_shape_priors: Any,
    exponent_by_channel_class: Any,
    modality_weights: Any,
    prior_factor: float,
    prior_coefficient: float,
    smoothing: float,
    prior_scale: float,
    intensity_field: Any | None = None,
    psf_bias: Any | None = None,
) -> ToadsInitialStatisticsMembershipContinuation:
    """Refresh stats on the P18cyn thinned labels, then run exact memberships.

    The first membership map is taken from the startup composition. Refreshed
    deviations feed the next membership stage as class variance weights and
    also drive source centroid-similarity coefficients. The image variance,
    modality weights, priors, exponents, transformed shapes, and relation
    labels are explicit fixed inputs, as in the selected Java call boundary.
    """
    if not isinstance(startup_and_thinning, ToadsInitialThinningComposition):
        raise TypeError("startup_and_thinning must be a selected startup composition")
    images = np.asarray(image_values, dtype=np.float32)
    labels = startup_and_thinning.thinning.segmentation
    memberships = startup_and_thinning.membership_pass.memberships
    if images.ndim != 4 or images.shape[:3] != labels.shape:
        raise ValueError("image_values must match the thinned segmentation grid")
    if memberships.shape != labels.shape + (11,):
        raise ValueError("startup memberships must have the selected 11-class layout")
    field = None if intensity_field is None else np.asarray(intensity_field, dtype=np.float32)
    if field is not None and field.shape != images.shape:
        raise ValueError("intensity_field must match image_values")

    centers = toads_default_exact_centroid_update(
        images,
        memberships,
        initial_centroids,
        category_exponents,
        intensity_field=field,
    )
    deviations = toads_default_exact_class_variances(
        images,
        memberships,
        centers,
        labels,
        intensity_field=field,
    )
    similarity = toads_optimized_centroid_similarity(
        centers, deviations, image_variance_weights, prior_scale
    )
    next_membership = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=memberships,
        centroids=centers,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        class_variance_weights=deviations,
        image_variances=image_variance_weights,
        transformed_shapes=transformed_shapes,
        closest_labels=startup_and_thinning.startup.relation_state.closest_labels,
        smoothing_priors=smoothing_priors,
        centroid_similarity=similarity,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        intensity_field=field,
        psf_bias=psf_bias,
    )
    return ToadsInitialStatisticsMembershipContinuation(
        startup_and_thinning=startup_and_thinning,
        centroids=centers,
        class_standard_deviations=deviations,
        centroid_similarity=similarity,
        membership_pass=next_membership,
        stage_order=(
            *startup_and_thinning.stage_order,
            "computeExactCentroids",
            "computeExactVariances",
            "computeCentroidSimilarity",
            "computeExactMemberships",
        ),
    )
