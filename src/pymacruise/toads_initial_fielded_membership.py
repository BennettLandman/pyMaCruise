"""Fixed-state composition of TOADS field attachment and initial membership.

The selected default 11-class wrapper attaches each correction field before
its startup ``computeExactMemberships`` call. This adapter connects those
per-channel references to the already validated 11-class membership kernel.
Upstream relation, boundary, closest-label, variance, image-variance, shape,
and similarity calculations remain explicit fixed inputs.
"""

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
from .toads_inhomogeneity_state import ToadsInhomogeneityFieldState


@dataclass(frozen=True)
class ToadsInitialFieldedMembership:
    """Initial membership output plus the exact attached-field references."""

    field_state: ToadsInhomogeneityFieldState
    membership_pass: ToadsDefaultAtlasMembershipPass
    stage_order: tuple[str, ...]


def toads_initial_default_atlas_membership_after_field_attachment(
    image_values: Any,
    initial_memberships: Any,
    centroids: Any,
    exponent_by_channel_class: Any,
    modality_weights: Any,
    class_variance_weights: Any,
    image_variances: Any,
    transformed_shapes: Any,
    closest_labels: Any,
    smoothing_priors: Any,
    centroid_similarity: Any,
    atlas_shape_priors: Any,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    *,
    channel_fields: Sequence[Any],
    class_types: Sequence[str] = DEFAULT_ATLAS_CLASS_TYPES,
    psf_bias: Any | None = None,
) -> ToadsInitialFieldedMembership:
    """Attach channel fields and run the selected fixed-state 11-class pass.

    Arrays are supplied in true RAS voxel order. ``channel_fields`` contains
    one fitted ``(x,y,z)`` field per image channel. Each field is assigned by
    reference, matching the selected Java setter, then values are assembled
    as ``(x,y,z,channel)`` for the existing Java-matched membership kernel.

    The source wrapper has already called ``computeRelations``,
    ``computeBoundaryDistances``, ``computeClosestLabels``, and
    ``computeVariances`` before its first ``computeMemberships`` call. Their
    outputs, along with transformed atlas shapes and statistics, are fixed
    inputs here. This function does not implement those upstream stages.
    """
    images = np.asarray(image_values, dtype=np.float32)
    if images.ndim != 4 or images.shape[-1] < 1:
        raise ValueError("image_values must have shape (x,y,z,channel)")
    if len(channel_fields) != images.shape[-1]:
        raise ValueError("channel_fields must contain one field per image channel")

    state = ToadsInhomogeneityFieldState(channel_count=images.shape[-1])
    for channel, field in enumerate(channel_fields):
        values = np.asarray(field, dtype=np.float32)
        if values.shape != images.shape[:3]:
            raise ValueError("each channel field must match the image spatial grid")
        state.add_inhomogeneity_correction(field, channel)
    if not state.use_field:
        raise ValueError("selected fielded membership requires an attached field")

    # Build the channel-last numerical array from the stored references. This
    # preserves source values while keeping the field-state setter semantics
    # separately observable on the returned object.
    field_volume: NDArray[np.float32] = np.stack(
        [np.asarray(value, dtype=np.float32) for value in state.field_by_channel],
        axis=-1,
    )
    membership = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=initial_memberships,
        centroids=centroids,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        image_variances=image_variances,
        transformed_shapes=transformed_shapes,
        closest_labels=closest_labels,
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
    return ToadsInitialFieldedMembership(
        field_state=state,
        membership_pass=membership,
        stage_order=tuple(
            [f"addInhomogeneityCorrection[channel={channel}]" for channel in range(images.shape[-1])]
            + [
            "fixed:computeRelations",
            "fixed:computeBoundaryDistances",
            "fixed:computeClosestLabels",
            "fixed:computeVariances",
            "computeExactMemberships",
            ]
        ),
    )
