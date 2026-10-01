"""Compose selected base-TOADS startup preparation and initial membership.

This adapter wires only the selected startup order: attach each channel's
field reference, initialize relation/boundary/closest-label state from the
constructor's zero boundary map, then pass the fixed statistics and
class-aligned transformed shapes into the 11-class exact-membership kernel.
It does not estimate those statistics, register/sample the atlas, or iterate
the TOADS segmentation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

from .toads_default_atlas_membership import (
    DEFAULT_ATLAS_CLASS_TYPES,
    ToadsDefaultAtlasMembershipPass,
    toads_default_atlas_membership_pass,
)
from .toads_inhomogeneity_state import ToadsInhomogeneityFieldState
from .toads_initial_relation_state import toads_initialize_relation_state
from .toads_relations import ToadsRelationState


@dataclass(frozen=True)
class ToadsPreparedInitialMembership:
    """Selected startup state and fixed-input membership output."""

    field_state: ToadsInhomogeneityFieldState
    relation_state: ToadsRelationState
    membership_pass: ToadsDefaultAtlasMembershipPass
    stage_order: tuple[str, ...]


def toads_prepared_initial_default_atlas_membership(
    image_values: Any,
    segmentation: Any,
    initial_memberships: Any,
    centroids: Any,
    exponent_by_channel_class: Any,
    modality_weights: Any,
    class_variance_weights: Any,
    image_variances: Any,
    transformed_shapes: Any,
    smoothing_priors: Any,
    centroid_similarity: Any,
    atlas_shape_priors: Any,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    *,
    channel_fields: Sequence[Any],
    limit: int,
    class_types: Sequence[str] = DEFAULT_ATLAS_CLASS_TYPES,
    psf_bias: Any | None = None,
) -> ToadsPreparedInitialMembership:
    """Run the default 11-class startup composition in selected source order.

    ``image_values`` is `(x,y,z,channel)` and all scalar fields, labels,
    memberships, and transformed shape maps use the same RAS voxel grid.
    `transformed_shapes[..., class]` follows the fixed default atlas descriptor
    order. Statistics, centroid similarity, shape samples and priors are
    explicit fixed inputs; `computeVariances`, image-variance estimation,
    similarity calculation, and atlas shape sampling are not recomputed here.

    Each field is assigned once through the selected Java setter semantics,
    relation preparation follows in the source's three-method order, and the
    resulting closest-label map is passed directly to the existing exact
    membership kernel. This avoids calling the fielded-membership adapter,
    which also performs attachment and would therefore repeat/reorder it.
    """
    images = np.asarray(image_values, dtype=np.float32)
    labels = np.asarray(segmentation)
    memberships = np.asarray(initial_memberships, dtype=np.float32)
    if images.ndim != 4 or images.shape[-1] < 1:
        raise ValueError("image_values must have shape (x,y,z,channel)")
    if labels.shape != images.shape[:3] or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer label grid matching image_values")
    if memberships.ndim != 4 or memberships.shape[:3] != images.shape[:3]:
        raise ValueError("initial_memberships must match the image grid and have a class axis")
    if memberships.shape[-1] != len(DEFAULT_ATLAS_CLASS_TYPES):
        raise ValueError("initial_memberships must use the selected default 11-class order")
    if len(channel_fields) != images.shape[-1]:
        raise ValueError("channel_fields must contain one field per image channel")
    if len(class_types) != len(DEFAULT_ATLAS_CLASS_TYPES) or tuple(class_types) != DEFAULT_ATLAS_CLASS_TYPES:
        raise ValueError("class_types must match the selected default 11-class atlas order")

    if any(np.asarray(field).shape != images.shape[:3] for field in channel_fields):
        raise ValueError("each channel field must match the image spatial grid")

    # AlgorithmMultipleToads attaches fields before initializing relation state.
    field_state = ToadsInhomogeneityFieldState(channel_count=images.shape[-1])
    stage_order: list[str] = []
    for channel, field in enumerate(channel_fields):
        field_state.add_inhomogeneity_correction(field, channel)
        stage_order.append(f"addInhomogeneityCorrection[channel={channel}]")

    relation_state = toads_initialize_relation_state(labels, memberships, limit)
    stage_order.extend(
        ("computeRelations", "computeBoundaryDistances", "computeClosestLabels")
    )

    # The source invokes computeVariances before computeExactMemberships. The
    # fixed-input adapter receives these results and the per-voxel atlas shapes.
    stage_order.append("fixed:computeVariances")
    field_volume = np.stack(
        [np.asarray(value, dtype=np.float32) for value in field_state.field_by_channel],
        axis=-1,
    )
    membership_pass = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=memberships,
        centroids=centroids,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        image_variances=image_variances,
        transformed_shapes=transformed_shapes,
        closest_labels=relation_state.closest_labels,
        smoothing_priors=smoothing_priors,
        centroid_similarity=centroid_similarity,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        class_types=class_types,
        intensity_field=field_volume,
        psf_bias=psf_bias,
    )
    stage_order.extend(
        (
            "fixed:computeImageVariance+centroidSimilarity+transformedShapes",
            "computeExactMemberships",
        )
    )
    return ToadsPreparedInitialMembership(
        field_state=field_state,
        relation_state=relation_state,
        membership_pass=membership_pass,
        stage_order=tuple(stage_order),
    )
