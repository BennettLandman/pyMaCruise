"""Compose TOADS initial statistics/membership with the first thinning pass."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_initial_relation_state import toads_initialize_relation_state
from .toads_initial_statistics_membership import (
    ToadsInitialStatisticsMembership,
    compute_toads_initial_statistics_membership,
)
from .toads_relations import ToadsRelationState
from .toads_thinning import ToadsThinningPass, toads_propagate_thinning


TOADS_INITIAL_STATISTICS_THINNING_ORDER = (
    "computeRelations",
    "computeBoundaryDistances",
    "computeClosestLabels",
    "computeVariances",
    "computeImageVariance",
    "computeExactMemberships.centroidSimilarity",
    "computeExactMemberships.voxelSweep",
    "initThinningLabels",
    "propagateThinning",
)


@dataclass(frozen=True)
class ToadsInitialStatisticsThinning:
    """Initial TOADS membership result followed by the first thinning state."""

    relation_state: ToadsRelationState
    statistics_membership: ToadsInitialStatisticsMembership
    thinning: ToadsThinningPass
    stage_order: tuple[str, ...]


def compute_toads_initial_statistics_thinning(
    *,
    image_values: Any,
    initial_memberships: Any,
    centroids: Any,
    segmentation: Any,
    class_types: Sequence[str],
    exponent_by_channel_class: Any,
    modality_weights: Any,
    transformed_shapes: Any,
    smoothing_priors: Any,
    atlas_shape_priors: Any,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    relation_limit: int = 3,
    volume_ratio: float = 0.5,
    thinning_limit: float = 3.0,
    prior_scale: float = 0.05,
    intensity_field: Any | None = None,
    psf_bias: Any | None = None,
) -> ToadsInitialStatisticsThinning:
    """Compute relations, initial membership, and the selected first thinning.

    Atlas registration, shape sampling, and field fitting are fixed upstream
    inputs. All volume arrays share RAS ``(x,y,z)`` coordinates and descriptor
    class order. The source-compatible relation state supplies both the
    closest-label map for membership and class adjacency for thinning.
    """
    labels = np.asarray(segmentation)
    memberships = np.asarray(initial_memberships, dtype=np.float32)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    if labels.ndim != 3 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be a 3-D integer class-index grid")
    if memberships.ndim != 4 or memberships.shape[:3] != labels.shape:
        raise ValueError("initial_memberships must match segmentation and contain a class axis")
    if shapes.shape != memberships.shape or not np.isfinite(shapes).all():
        raise ValueError("transformed_shapes must be finite and match initial_memberships")
    relation_state = toads_initialize_relation_state(labels, memberships, relation_limit)
    statistics_membership = compute_toads_initial_statistics_membership(
        image_values=image_values,
        initial_memberships=memberships,
        centroids=centroids,
        segmentation=labels,
        class_types=class_types,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        transformed_shapes=shapes,
        closest_labels=relation_state.closest_labels,
        smoothing_priors=smoothing_priors,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        prior_scale=prior_scale,
        intensity_field=intensity_field,
        psf_bias=psf_bias,
    )
    thinning = toads_propagate_thinning(
        labels,
        statistics_membership.membership_pass.memberships,
        relation_state.relation.class_adjacency,
        shapes,
        volume_ratio=volume_ratio,
        limit=thinning_limit,
    )
    return ToadsInitialStatisticsThinning(
        relation_state,
        statistics_membership,
        thinning,
        TOADS_INITIAL_STATISTICS_THINNING_ORDER,
    )
