"""Compose one selected default-TOADS startup-to-membership state transition.

This wires the already source-matched exact initializer/thinning, centroid,
variance, similarity, and 11-class membership kernels. It treats image-level
variance weights, transformed shapes, closest relation labels, and descriptor
weights as fixed inputs; it is not the complete wrapper workflow.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_default_atlas_membership import (
    ToadsDefaultAtlasMembershipPass,
    toads_default_atlas_membership_pass,
)
from .toads_default_exact_centroids import toads_default_exact_centroid_update
from .toads_default_exact_variances import toads_default_exact_class_variances
from .toads_first_thinning import ToadsFirstThinningPass, toads_first_exact_thinning_pass
from .toads_math import toads_optimized_centroid_similarity


@dataclass(frozen=True)
class ToadsDefaultExactIterationSlice:
    """Fixed-input outputs at each stage of one selected exact iteration slice."""

    first_thinning: ToadsFirstThinningPass
    refreshed_centroids: NDArray[np.float32]
    class_standard_deviations: NDArray[np.float32]
    centroid_similarity: NDArray[np.float32]
    membership_pass: ToadsDefaultAtlasMembershipPass


def toads_default_exact_iteration_slice(
    image_values: ArrayLike,
    initial_centroids: ArrayLike,
    initial_segmentation: ArrayLike,
    initial_thinning_relations: ArrayLike,
    category_exponents: ArrayLike,
    modality_weights: ArrayLike,
    class_variance_weights: ArrayLike,
    image_variance_weights: ArrayLike,
    transformed_shapes: ArrayLike,
    closest_labels: ArrayLike,
    smoothing_priors: ArrayLike,
    atlas_shape_priors: ArrayLike,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    *,
    intensity_field: ArrayLike | None = None,
    crop_border: int = 10,
    volume_ratio: float = 0.5,
    thinning_limit: float = 3.0,
    centroid_smoothness: float = 0.1,
    prior_scale: float = 0.05,
) -> ToadsDefaultExactIterationSlice:
    """Run first exact thinning through one exact 11-class membership sweep.

    ``image_variance_weights`` are a fixed post-``computeImageVariance`` input
    and are used both for centroid similarity and the data-energy pass. The
    stage deliberately excludes image-wide variance estimation, field fitting,
    registration, further loop/stop/convergence logic, growing and final
    classification. ``transformed_shapes`` and ``closest_labels`` are also
    fixed RAS-grid inputs describing the shape/relation state for this pass.
    """
    images = np.asarray(image_values, dtype=np.float32)
    centers = np.asarray(initial_centroids, dtype=np.float32)
    labels = np.asarray(initial_segmentation)
    field = None if intensity_field is None else np.asarray(intensity_field, dtype=np.float32)
    if images.ndim != 4 or images.shape[-1] != 2 or not np.isfinite(images).all():
        raise ValueError("image_values must be finite with shape (x,y,z,2)")
    if centers.shape != (2, 11) or not np.isfinite(centers).all():
        raise ValueError("initial_centroids must be finite with shape (2,11)")
    if labels.shape != images.shape[:3] or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("initial_segmentation must be an integer map matching image_values")
    if field is not None and (field.shape != images.shape or not np.isfinite(field).all()):
        raise ValueError("intensity_field must be finite and match image_values")

    fixed_shapes = np.asarray(transformed_shapes, dtype=np.float32)
    fixed_labels = np.asarray(closest_labels)
    if fixed_shapes.shape != images.shape[:3] + (11,) or not np.isfinite(fixed_shapes).all():
        raise ValueError("transformed_shapes must be finite and match the 11-class image grid")
    if fixed_labels.shape != images.shape[:3] + (4,) or not np.issubdtype(fixed_labels.dtype, np.integer):
        raise ValueError("closest_labels must be an integer four-slot map matching image_values")

    first = toads_first_exact_thinning_pass(
        np.moveaxis(images, -1, 0),
        centers,
        labels,
        initial_thinning_relations,
        np.zeros(images.shape[:3] + (11,), dtype=np.float32),
        crop_border=crop_border,
        volume_ratio=volume_ratio,
        limit=thinning_limit,
    )
    centroids = toads_default_exact_centroid_update(
        images,
        first.memberships,
        centers,
        category_exponents,
        intensity_field=field,
        smoothness=centroid_smoothness,
    )
    deviations = toads_default_exact_class_variances(
        images,
        first.memberships,
        centroids,
        first.thinning.segmentation,
        intensity_field=field,
    )
    image_var = np.asarray(image_variance_weights, dtype=np.float32)
    similarity = toads_optimized_centroid_similarity(
        centroids, deviations, image_var, prior_scale=prior_scale
    )
    membership = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=first.memberships,
        centroids=centroids,
        exponent_by_channel_class=_expand_default_exponents(category_exponents),
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        image_variances=image_var,
        transformed_shapes=fixed_shapes,
        closest_labels=fixed_labels,
        smoothing_priors=smoothing_priors,
        centroid_similarity=similarity,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        intensity_field=field,
    )
    return ToadsDefaultExactIterationSlice(first, centroids, deviations, similarity, membership)


def _expand_default_exponents(values: Any) -> NDArray[np.float32]:
    category = np.asarray(values, dtype=np.float32)
    if category.shape != (2, 3) or not np.isfinite(category).all():
        raise ValueError("category_exponents must be finite with shape (2,CSF/GM/WM)")
    class_category = (0, 0, 1, 1, 1, 1, 1, 2, 2, 2)
    expanded = np.empty((2, 11), dtype=np.float32)
    expanded[:, 0] = np.float32(1.0)
    for class_index, category_index in enumerate(class_category, start=1):
        expanded[:, class_index] = category[:, category_index]
    return expanded
