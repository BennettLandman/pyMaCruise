"""Source-faithful TOADS initialization with unity correction fields."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_initial_statistics_membership import ToadsInitialStatisticsMembership
from .toads_registered_initial_membership import compute_registered_toads_initial_membership
from .toads_thinning import ToadsThinningPass, toads_propagate_thinning


@dataclass(frozen=True)
class ToadsRegisteredUnityFieldStartup:
    """Selected registered startup through its initial membership and thinning."""

    fields: NDArray[np.float32]
    degree_state: int
    energy_ratio: np.float32
    statistics_membership: ToadsInitialStatisticsMembership
    thinning: ToadsThinningPass
    stage_order: tuple[str, ...]


def compute_registered_toads_unity_field_startup(
    alignment_run: Any,
    *,
    exponent_by_channel_class: Any,
    modality_weights: Any,
    smoothing_priors: Any,
    atlas_shape_priors: Any,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    prior_scale: float = 0.05,
    volume_ratio: float = 0.5,
    thinning_limit: float = 3.0,
    psf_bias: Any | None = None,
    record_thinning_items: bool = True,
) -> ToadsRegisteredUnityFieldStartup:
    """Run the selected initial pass without prematurely fitting a field.

    ``AlgorithmMultipleToads`` constructs each correction field as unity and
    starts with ``deg=0`` and ``Eratio=1``. It computes initial memberships
    before entering the loop. The first loop stage is thinning; its degree-zero
    field update is a no-op. This composition returns the initial membership
    pass and that first thinning result with those carried control values.
    Atlas and segmentation coefficients remain explicit inputs.
    """
    try:
        recropped = np.asarray(alignment_run.template_refresh.recropped_images, dtype=np.float32)
        segmentation = np.asarray(alignment_run.segmentation_constructor.segmentation)
        initial_memberships = np.asarray(
            alignment_run.segmentation_initialization.memberships, dtype=np.float32
        )
        shapes = np.asarray(alignment_run.transformed_shape_priors, dtype=np.float32)
        adjacency = np.asarray(alignment_run.relation_state.relation.class_adjacency)
    except AttributeError as exc:
        raise TypeError("alignment_run is missing selected registered TOADS state") from exc
    if recropped.ndim != 4 or recropped.shape[0] < 1:
        raise ValueError("recropped_images must have shape (channel,x,y,z)")
    images = np.moveaxis(recropped, 0, -1)
    if segmentation.shape != images.shape[:3]:
        raise ValueError("constructor segmentation must match the recropped image grid")
    if initial_memberships.shape != images.shape[:3] + (11,):
        raise ValueError("initial memberships must match the selected 11-class image grid")
    if shapes.shape != initial_memberships.shape:
        raise ValueError("transformed shape priors must match initial memberships")
    if adjacency.shape != (11, 11):
        raise ValueError("registered class adjacency must match the selected 11 classes")

    fields = np.ones(images.shape, dtype=np.float32)
    membership = compute_registered_toads_initial_membership(
        alignment_run,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        smoothing_priors=smoothing_priors,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        prior_scale=prior_scale,
        intensity_field=fields,
        psf_bias=psf_bias,
    )
    thinning = toads_propagate_thinning(
        segmentation,
        membership.membership_pass.memberships,
        adjacency,
        shapes,
        volume_ratio=volume_ratio,
        limit=thinning_limit,
        record_processed_items=record_thinning_items,
    )
    return ToadsRegisteredUnityFieldStartup(
        fields=fields,
        degree_state=0,
        energy_ratio=np.float32(1.0),
        statistics_membership=membership,
        thinning=thinning,
        stage_order=(
            "registered_startup.initializeCorrectionFieldsToUnity",
            "registered_startup.computeVariances",
            "registered_startup.computeImageVariance",
            "registered_startup.computeCentroidSimilarity",
            "registered_startup.computeExactMemberships",
            "iteration1.propagateThinning",
        ),
    )
