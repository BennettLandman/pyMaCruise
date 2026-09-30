"""Compose the next fixed-statistics TOADS thinning and membership stages.

This adapter consumes the continuing post-registration result: its grown
segmentation, refreshed class relations/closest labels, current rigid
transform, and current membership maps. It samples the eleven atlas priors at
that transform, runs the selected source thinning queue, and immediately
runs the next exact-membership sweep with explicitly supplied statistics.
It does not estimate centroids/variances, register the atlas, or decide
whether another outer iteration should run. Paper context: Bazin and Pham
(2007, 2008), ``docs/papers/toads.md``. This connects one selected continuing
transition to the next tissue-class membership sweep; it is not convergence
or whole-algorithm validation.
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
from .toads_post_registration import ToadsPostRegistrationBranch
from .toads_shapes import sample_toad_shape_prior
from .toads_thinning import ToadsThinningPass, toads_propagate_thinning


@dataclass(frozen=True)
class ToadsNextDefaultMembershipPass:
    """Outputs and explicit order for one continuing-loop composition."""

    transformed_shapes: NDArray[np.float32]
    thinning: ToadsThinningPass
    membership_pass: ToadsDefaultAtlasMembershipPass
    stage_order: tuple[str, ...]


def _sample_current_transform_shapes(
    atlas_shapes: NDArray[np.float32],
    branch: ToadsPostRegistrationBranch,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
) -> NDArray[np.float32]:
    handoff = branch.handoff
    spatial = handoff.membership_pass.memberships.shape[:3]
    if len(handoff.shape_bounds) != 11 or atlas_shapes.shape[-1] != 11:
        raise ValueError("the continuing default-atlas handoff must provide eleven shape bounds and classes")
    transform = tuple(float(value) for value in handoff.transform)
    result = np.empty(spatial + (11,), dtype=np.float32)
    for class_index, bounds in enumerate(handoff.shape_bounds):
        template_shape = atlas_shapes[..., class_index]
        for x in range(spatial[0]):
            for y in range(spatial[1]):
                for z in range(spatial[2]):
                    result[x, y, z, class_index] = sample_toad_shape_prior(
                        template_shape,
                        (x, y, z),
                        image_dimensions_xyz=spatial,
                        image_resolutions_mm=image_resolutions_mm,
                        template_resolutions_mm=template_resolutions_mm,
                        transform=transform,
                        class_index=class_index,
                        bbox_low_xyz=bounds.lower_xyz,
                        bbox_high_xyz=bounds.upper_xyz,
                    ).value
    return result


def toads_next_default_membership_pass(
    branch: ToadsPostRegistrationBranch,
    *,
    image_values: Any,
    centroids: Any,
    exponent_by_channel_class: Any,
    modality_weights: Any,
    class_variance_weights: Any,
    image_variances: Any,
    smoothing_priors: Any,
    centroid_similarity: Any,
    atlas_shape_priors: Any,
    atlas_shapes: Any,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    volume_ratio: float,
    thinning_limit: float,
    intensity_field: Any | None = None,
    psf_bias: Any | None = None,
) -> ToadsNextDefaultMembershipPass:
    """Run ``propagateThinning`` then the next selected default membership pass.

    The branch must be the non-stopping P18czi result. Its refreshed relation
    adjacency and closest-label map are consumed directly. Class centroids,
    deviations/variance weights, image variances, class exponents, and all
    other statistics are fixed call inputs, preserving the wrapper stage
    boundary while excluding statistical estimation.
    """
    if not isinstance(branch, ToadsPostRegistrationBranch):
        raise TypeError("branch must be a ToadsPostRegistrationBranch")
    if branch.control.stop or branch.growing is None or branch.relations is None:
        raise ValueError("a continuing post-registration branch with growth and refreshed relations is required")
    if branch.control.update_classification:
        raise ValueError("a stopping/classification branch cannot continue into TOADS thinning")

    images = np.asarray(image_values, dtype=np.float32)
    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    memberships = branch.handoff.membership_pass.memberships
    if images.ndim != 4 or images.shape[:3] != memberships.shape[:3] or images.shape[-1] != 2:
        raise ValueError("image_values must be two-channel data matching the handoff subject grid")
    if shapes.ndim != 4 or shapes.shape[-1] != 11 or not np.isfinite(shapes).all():
        raise ValueError("atlas_shapes must be a finite (template_x, template_y, template_z, 11) array")
    if branch.relations.relation.class_adjacency.shape != (11, 11):
        raise ValueError("refreshed relation state must contain the default 11-class adjacency")
    if branch.relations.closest_labels.shape != memberships.shape[:3] + (4,):
        raise ValueError("refreshed closest labels must match the subject grid")

    transformed_shapes = _sample_current_transform_shapes(
        shapes, branch, image_resolutions_mm, template_resolutions_mm
    )
    thinned = toads_propagate_thinning(
        branch.growing.segmentation,
        memberships,
        branch.relations.relation.class_adjacency,
        transformed_shapes,
        volume_ratio,
        thinning_limit,
    )
    membership = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=memberships,
        centroids=centroids,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        image_variances=image_variances,
        transformed_shapes=transformed_shapes,
        closest_labels=branch.relations.closest_labels,
        smoothing_priors=smoothing_priors,
        centroid_similarity=centroid_similarity,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        intensity_field=intensity_field,
        psf_bias=psf_bias,
    )
    return ToadsNextDefaultMembershipPass(
        transformed_shapes,
        thinned,
        membership,
        (
            "continuing_post_registration_result",
            "propagateThinning",
            "computeExactMemberships",
        ),
    )
