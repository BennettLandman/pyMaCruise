"""Compose selected default TOADS startup state into its first thinning pass.

This bounded composition starts from fixed upstream segmentation,
memberships, and image statistics. It attaches the supplied correction fields,
refreshes relation/boundary/closest-label state, evaluates the exact
11-class membership pass, and sends that result into the selected first
thinning kernel. Atlas registration, estimation of the supplied statistics,
and later TOADS iterations remain outside this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_default_atlas_membership import ToadsDefaultAtlasMembershipPass
from .toads_initial_prepared_membership import (
    ToadsPreparedInitialMembership,
    toads_prepared_initial_default_atlas_membership,
)
from .toads_thinning import ToadsThinningPass, toads_propagate_thinning


@dataclass(frozen=True)
class ToadsInitialThinningComposition:
    """Observed startup states through the selected first thinning result."""

    startup: ToadsPreparedInitialMembership
    thinning: ToadsThinningPass
    stage_order: tuple[str, ...]

    @property
    def membership_pass(self) -> ToadsDefaultAtlasMembershipPass:
        return self.startup.membership_pass


def toads_initial_fielded_membership_to_thinning(
    image_values: ArrayLike,
    segmentation: ArrayLike,
    initial_memberships: ArrayLike,
    centroids: ArrayLike,
    exponent_by_channel_class: ArrayLike,
    modality_weights: ArrayLike,
    class_variance_weights: ArrayLike,
    image_variances: ArrayLike,
    transformed_shapes: ArrayLike,
    smoothing_priors: ArrayLike,
    centroid_similarity: ArrayLike,
    atlas_shape_priors: ArrayLike,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    *,
    channel_fields: Sequence[ArrayLike],
    limit: int,
    volume_ratio: float = 0.5,
    thinning_limit: float = 3.0,
    psf_bias: ArrayLike | None = None,
    initial_skeleton: ArrayLike | None = None,
) -> ToadsInitialThinningComposition:
    """Run startup relation/membership state and one thinning pass in order.

    All spatial arrays are on the same true-RAS voxel grid. Class statistics,
    image variance, centroid similarity, atlas priors, and transformed shape
    samples are supplied fixed inputs. In particular, this function does not
    create registration transforms or sample atlas shapes from a T1.
    """
    labels = np.asarray(segmentation)
    memberships = np.asarray(initial_memberships)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    if labels.ndim != 3 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be a 3-D integer label map")
    if memberships.ndim != 4 or memberships.shape[:3] != labels.shape:
        raise ValueError("initial_memberships must match segmentation and have a class axis")
    if shapes.shape != memberships.shape or not np.isfinite(shapes).all():
        raise ValueError("transformed_shapes must be finite and match the membership grid")

    startup = toads_prepared_initial_default_atlas_membership(
        image_values=image_values,
        segmentation=labels,
        initial_memberships=memberships,
        centroids=centroids,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        image_variances=image_variances,
        transformed_shapes=shapes,
        smoothing_priors=smoothing_priors,
        centroid_similarity=centroid_similarity,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        channel_fields=channel_fields,
        limit=limit,
        psf_bias=psf_bias,
    )
    thinning = toads_propagate_thinning(
        labels,
        startup.membership_pass.memberships,
        startup.relation_state.relation.class_adjacency,
        shapes,
        volume_ratio=volume_ratio,
        limit=thinning_limit,
        # The Java constructor starts with a zero-filled byte array. Replays
        # of harnesses with preinitialized state can provide that seed.
        initial_skeleton=initial_skeleton,
    )
    return ToadsInitialThinningComposition(
        startup=startup,
        thinning=thinning,
        stage_order=startup.stage_order + ("initThinningLabels", "propagateThinning"),
    )
