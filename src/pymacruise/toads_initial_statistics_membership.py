"""Compose source-ordered TOADS statistics and the initial exact pass."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from numpy.typing import NDArray

from .toads_default_atlas_membership import (
    DEFAULT_ATLAS_CLASS_TYPES,
    ToadsDefaultAtlasMembershipPass,
    toads_default_atlas_membership_pass,
)
from .toads_segmentation_initial_image_variance import (
    ToadsInitialImageVariance,
    initialize_toads_image_variance,
)
from .toads_segmentation_initial_similarity import (
    ToadsInitialCentroidSimilarity,
    initialize_toads_centroid_similarity,
)
from .toads_segmentation_initial_variance import (
    ToadsInitialClassVariances,
    initialize_toads_class_variances,
)


TOADS_INITIAL_STATISTICS_MEMBERSHIP_ORDER = (
    "segmentation.computeVariances",
    "segmentation.computeImageVariance",
    "segmentation.computeExactMemberships.centroidSimilarity",
    "segmentation.computeExactMemberships.voxelSweep",
)


@dataclass(frozen=True)
class ToadsInitialStatisticsMembership:
    """Initial TOADS statistics and exact membership output for fixed priors."""

    class_variances: ToadsInitialClassVariances
    image_variance: ToadsInitialImageVariance
    centroid_similarity: ToadsInitialCentroidSimilarity
    membership_pass: ToadsDefaultAtlasMembershipPass
    stage_order: tuple[str, ...]


def compute_toads_initial_statistics_membership(
    *,
    image_values: Any,
    initial_memberships: Any,
    centroids: Any,
    segmentation: Any,
    class_types: Sequence[str],
    class_variance_weights: Any,
    exponent_by_channel_class: Any,
    modality_weights: Any,
    transformed_shapes: Any,
    closest_labels: Any,
    smoothing_priors: Any,
    atlas_shape_priors: Any,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    prior_scale: float = 0.05,
    intensity_field: Any | None = None,
    psf_bias: Any | None = None,
) -> ToadsInitialStatisticsMembership:
    """Refresh selected base-TOADS statistics, then run one exact pass.

    The input arrays represent the source state immediately before
    ``computeVariances``. Refreshed standard deviations feed centroid
    similarity; fixed atlas Intensity Variance Factors feed memberships.
    Shape samples, closest labels, priors, and PowerTable exponents are fixed
    upstream inputs, so this adapter composes only the selected
    variance/similarity/membership sequence. Arrays are in RAS
    ``(x,y,z[,channel])`` order and descriptor class order.
    """
    types = tuple(str(value).upper() for value in class_types)
    if types != DEFAULT_ATLAS_CLASS_TYPES:
        raise ValueError("class_types must match the selected default 11-class atlas order")
    images = np.asarray(image_values, dtype=np.float32)
    memberships = np.asarray(initial_memberships, dtype=np.float32)
    labels = np.asarray(segmentation)
    centers = np.asarray(centroids, dtype=np.float32)
    if images.ndim != 4 or images.shape[-1] < 1 or min(images.shape[:3]) < 3:
        raise ValueError("image_values must have shape (x>=3,y>=3,z>=3,channel>=1)")
    if memberships.shape != images.shape[:3] + (len(types),):
        raise ValueError("initial_memberships must match the image grid and 11-class order")
    if labels.shape != images.shape[:3] or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer class-index grid matching image_values")
    if centers.shape != (images.shape[-1], len(types)):
        raise ValueError("centroids must have shape (channel,11)")
    if not np.isfinite(images).all() or not np.isfinite(memberships).all() or not np.isfinite(centers).all():
        raise ValueError("images, memberships, and centroids must be finite")

    class_variances = initialize_toads_class_variances(
        images_xyz_channels=images,
        memberships=memberships,
        centroids=centers,
        segmentation=labels,
        class_types=types,
    )
    image_variance = initialize_toads_image_variance(
        images_xyz_channels=images,
        memberships=memberships,
        centroids=centers,
        class_types=types,
    )
    similarity = initialize_toads_centroid_similarity(
        centroids=centers,
        class_standard_deviations=class_variances.standard_deviations,
        image_variance_weights=image_variance.image_variance_weights,
        class_types=types,
        prior_scale=prior_scale,
    )
    membership_pass = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=memberships,
        centroids=centers,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        image_variances=image_variance.image_variance_weights,
        transformed_shapes=transformed_shapes,
        closest_labels=closest_labels,
        smoothing_priors=smoothing_priors,
        centroid_similarity=similarity.similarity,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        class_types=types,
        intensity_field=intensity_field,
        psf_bias=psf_bias,
    )
    return ToadsInitialStatisticsMembership(
        class_variances,
        image_variance,
        similarity,
        membership_pass,
        TOADS_INITIAL_STATISTICS_MEMBERSHIP_ORDER,
    )
