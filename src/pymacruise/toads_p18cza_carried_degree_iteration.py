"""Bounded next-loop degree-one continuation with source-carried degree state.

The selected ``AlgorithmMultipleToads`` loop keeps its polynomial degree
between iterations. This adapter composes one P18cyz thinning output through
the already validated degree-one field, statistics, and membership kernels,
while explicitly applying the prior degree instead of restarting the ramp.
It is a single-channel synthetic-state slice, not a complete TOADS workflow.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .toads_default_atlas_membership import (
    ToadsDefaultAtlasMembershipPass,
    toads_default_atlas_membership_pass,
)
from .toads_default_exact_centroids import toads_default_exact_centroid_update
from .toads_default_exact_variances import toads_default_exact_class_variances
from .toads_inhomogeneity_exact_fit import (
    fit_toads_exact_chebyshev_field_degree1,
)
from .toads_inhomogeneity_transition import toads_next_inhomogeneity_degree
from .toads_math import toads_optimized_centroid_similarity, toads_optimized_image_variance
from .toads_post_registration_next_thinning import ToadsPostRegistrationNextThinning
from .toads_post_thinning_exact_field import ToadsPostThinningExactField
from .toads_post_thinning_exact_iteration import ToadsPostThinningExactIteration


@dataclass(frozen=True)
class ToadsP18czaCarriedDegreeIteration:
    """One exact iteration plus the carried degree and input branch state."""

    next_thinning: ToadsPostRegistrationNextThinning
    previous_degree: int
    degree: int
    energy_ratio: np.float32
    exact_iteration: ToadsPostThinningExactIteration
    stage_order: tuple[str, ...]


def continue_p18cyz_thinning_through_carried_degree1_iteration(
    next_thinning: ToadsPostRegistrationNextThinning,
    image_values: Any,
    initial_centroids: Any,
    field: Any,
    *,
    previous_degree: int,
    polynomial_degree: int,
    low_energy_ratio: float,
    category_exponents: Any,
    exponent_by_channel_class: Any,
    modality_weights: Any,
    class_variance_weights: Any,
    transformed_shapes: Any,
    closest_labels: Any,
    smoothing_priors: Any,
    atlas_shape_priors: Any,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    prior_scale: float,
    centroid_smoothness: float = 0.1,
) -> ToadsP18czaCarriedDegreeIteration:
    """Apply source loop order to a continuing branch with explicit degree state.

    The only supported scheduled degree here is one. The energy ratio is the
    just-completed P18cyx membership/control ratio consumed after the P18cyz
    thinning pass. ``previous_degree`` must reflect the earlier loop; the
    source increments it by at most one when ``ratio < low_energy_ratio``.
    """
    if not isinstance(next_thinning, ToadsPostRegistrationNextThinning):
        raise TypeError("next_thinning must be a selected post-registration thinning result")
    branch = next_thinning.branch
    if branch.control.stop or branch.growing is None or branch.relations is None:
        raise ValueError("carried-degree iteration requires a continuing registered branch")
    ratio = np.float32(branch.control.energy_ratio)
    degree = toads_next_inhomogeneity_degree(
        previous_degree, polynomial_degree, ratio, low_energy_ratio
    )
    if degree != 1:
        raise ValueError("this bounded continuation supports only source-scheduled degree one")

    image = np.asarray(image_values, dtype=np.float32)
    memberships = np.asarray(branch.handoff.membership_pass.memberships, dtype=np.float32)
    labels = np.asarray(next_thinning.thinning.segmentation)
    old_centers = np.asarray(initial_centroids, dtype=np.float32)
    category = np.asarray(category_exponents, dtype=np.float32)
    if image.ndim != 3 or image.shape != memberships.shape[:3]:
        raise ValueError("image_values must be one 3-D channel on the membership grid")
    if old_centers.shape != (11,):
        raise ValueError("initial_centroids must contain the selected 11-class values")
    if category.shape != (3,):
        raise ValueError("category_exponents must contain CSF, GM, and WM values")

    field_fit = fit_toads_exact_chebyshev_field_degree1(
        image, memberships, labels, old_centers, field
    )
    field_update = ToadsPostThinningExactField(
        next_thinning=next_thinning,
        degree=degree,
        fit=field_fit,
        stage_order=(
            *next_thinning.stage_order,
            "next_iteration.computeExact3DPolynomialField.degree1",
        ),
    )
    images = image[..., None]
    fields = field_fit.field[..., None]
    exponent_map = np.empty((1, 11), dtype=np.float32)
    exponent_map[0, 0] = np.float32(1.0)
    exponent_map[0, 1:3] = category[0]
    exponent_map[0, 3:8] = category[1]
    exponent_map[0, 8:11] = category[2]
    supplied_exponents = np.asarray(exponent_by_channel_class, dtype=np.float32)
    if supplied_exponents.shape != (1, 11) or not np.array_equal(supplied_exponents, exponent_map):
        raise ValueError("exponent_by_channel_class must match the one-channel default 11-class profile")

    centers = toads_default_exact_centroid_update(
        images,
        memberships,
        old_centers[None, :],
        category[None, :],
        intensity_field=fields,
        smoothness=centroid_smoothness,
    )
    deviations = toads_default_exact_class_variances(
        images, memberships, centers, labels, intensity_field=fields
    )
    image_variance_weights, image_std = toads_optimized_image_variance(
        images,
        memberships[..., 1:],
        centers[:, 1:],
        intensity_field=fields,
    )
    similarity = toads_optimized_centroid_similarity(
        centers, deviations, image_variance_weights, prior_scale
    )
    membership_pass = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=memberships,
        centroids=centers,
        exponent_by_channel_class=supplied_exponents,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        image_variances=image_variance_weights,
        transformed_shapes=transformed_shapes,
        closest_labels=closest_labels,
        smoothing_priors=smoothing_priors,
        centroid_similarity=similarity,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        intensity_field=fields,
    )
    exact = ToadsPostThinningExactIteration(
        field_update=field_update,
        centroids=centers,
        standard_deviations=deviations,
        image_variance_weights=image_variance_weights,
        image_standard_deviation=image_std,
        centroid_similarity=similarity,
        membership_pass=membership_pass,
        stage_order=(
            *field_update.stage_order,
            "next_iteration.computeExactCentroids",
            "next_iteration.computeExactVariances",
            "next_iteration.computeImageVariance",
            "next_iteration.computeCentroidSimilarity",
            "next_iteration.computeExactMemberships",
        ),
    )
    return ToadsP18czaCarriedDegreeIteration(
        next_thinning=next_thinning,
        previous_degree=int(previous_degree),
        degree=degree,
        energy_ratio=ratio,
        exact_iteration=exact,
        stage_order=exact.stage_order,
    )
