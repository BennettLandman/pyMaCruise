"""Selected single-channel TOADS post-thinning exact iteration continuation."""

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
from .toads_math import toads_optimized_centroid_similarity, toads_optimized_image_variance
from .toads_post_registration_next_thinning import ToadsPostRegistrationNextThinning
from .toads_post_thinning_exact_field import (
    ToadsPostThinningExactField,
    toads_post_thinning_exact_field_degree1,
)


@dataclass(frozen=True)
class ToadsPostThinningExactIteration:
    """Outputs of one selected degree-one field/statistics/membership sequence."""

    field_update: ToadsPostThinningExactField
    centroids: NDArray[np.float32]
    standard_deviations: NDArray[np.float32]
    image_variance_weights: NDArray[np.float32]
    image_standard_deviation: NDArray[np.float32]
    centroid_similarity: NDArray[np.float32]
    membership_pass: ToadsDefaultAtlasMembershipPass
    stage_order: tuple[str, ...]


def toads_post_thinning_exact_iteration_degree1(
    next_thinning: ToadsPostRegistrationNextThinning,
    image_values: Any,
    initial_centroids: Any,
    field: Any,
    *,
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
) -> ToadsPostThinningExactIteration:
    """Run selected next-loop field fit, statistics refresh, then exact memberships.

    This composes one channel and degree-one exact field fit only. The input
    T1 channel, previous membership maps, post-thinning labels, transform,
    shapes, closest-label map, and membership parameters are explicit fixed
    state. It does not run registration, estimate other inputs, or assemble
    the outer TOADS workflow. The field is updated in place as in the selected
    attached-array source operation.
    """
    if not isinstance(next_thinning, ToadsPostRegistrationNextThinning):
        raise TypeError("next_thinning must be a selected post-registration thinning result")
    if next_thinning.branch.control.stop:
        raise ValueError("an exact iteration continuation requires a continuing branch")
    image = np.asarray(image_values, dtype=np.float32)
    memberships = np.asarray(next_thinning.branch.handoff.membership_pass.memberships, dtype=np.float32)
    labels = np.asarray(next_thinning.thinning.segmentation)
    old_centers = np.asarray(initial_centroids, dtype=np.float32)
    category = np.asarray(category_exponents, dtype=np.float32)
    if image.ndim != 3 or image.shape != memberships.shape[:3]:
        raise ValueError("image_values must be one 3-D channel on the membership grid")
    if old_centers.shape != (11,):
        raise ValueError("initial_centroids must contain the selected 11-class values")
    if category.shape != (3,):
        raise ValueError("category_exponents must contain CSF, GM, and WM values")

    # AlgorithmMultipleToads.runAlgorithm orders the attached correction field
    # update before computeExactCentroids, computeExactVariances, and
    # computeMemberships. Keep the same mutable field object for each stage.
    field_update = toads_post_thinning_exact_field_degree1(
        next_thinning,
        image,
        old_centers,
        field,
    )
    images = image[..., None]
    fields = field_update.fit.field[..., None]
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
        images,
        memberships,
        centers,
        labels,
        intensity_field=fields,
    )
    # computeExactMemberships() refreshes computeImageVariance() before its
    # centroid-similarity and per-voxel membership calculations. MASK is not
    # an image-variance contributor in the selected switch; all other classes
    # in this descriptor are optimized tissue classes.
    image_variance_weights, image_std = toads_optimized_image_variance(
        images,
        memberships[..., 1:],
        centers[:, 1:],
        intensity_field=fields,
    )
    similarity = toads_optimized_centroid_similarity(
        centers,
        deviations,
        image_variance_weights,
        prior_scale,
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
    return ToadsPostThinningExactIteration(
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
