"""Compose the selected post-thinning LesionTOADS statistics refresh.

This adapter continues an already completed initial-thinning composition. It
applies the selected wrapper's ``computeCentroids`` → ``computeVariances`` →
``computeChannelWeights`` order using the existing float32 kernels. Inputs
such as the pre-thinning memberships and initialized black-hole centroid are
fixed state; this function does not infer or create missing classes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any
from collections.abc import Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_initial_thinning_composition import LesionInitialThinningComposition
from .lesion_statistics import (
    LesionToadsAuxiliaryCentroids,
    LesionToadsVariances,
    compute_lesion_toads_centroids,
    compute_lesion_toads_channel_weights,
    compute_lesion_toads_variances,
)
from .lesion_statistics_refresh import _source_direct_lesion_centroid


@dataclass(frozen=True)
class LesionPostThinningStatistics:
    """Statistics produced after the supplied first thinning pass."""

    centroids: LesionToadsAuxiliaryCentroids
    variances: LesionToadsVariances
    channel_weights: NDArray[np.float32]
    stage_order: tuple[str, ...]


def refresh_statistics_after_initial_thinning(
    composition: LesionInitialThinningComposition,
    images: ArrayLike,
    class_centroids: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    outlier_membership: ArrayLike,
    lesion_centroids: ArrayLike,
    black_hole_centroids: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    *,
    centroid_options: Mapping[str, Any] | None = None,
    channel_weight_options: Mapping[str, Any] | None = None,
    use_lesion_weight: bool = False,
    flair_channel: int = -1,
) -> LesionPostThinningStatistics:
    """Refresh statistics from the segmentation emitted by initial thinning.

    Membership arrays and initialized parameters remain fixed inputs, matching
    the wrapper at this call site. In particular, a synthetic profile without
    a BLACKHOLE class may still pass its initialized black-hole centroid and
    zero black-hole membership map; the helper does not synthesize that class.
    """
    if not isinstance(composition, LesionInitialThinningComposition):
        raise TypeError("composition must be a LesionInitialThinningComposition")
    if composition.stage_order[-2:] != (
        "propagateThinning:initThinningLabels", "propagateThinning:heapLoop"
    ):
        raise ValueError("composition must end at the selected initial thinning pass")
    centroid_kwargs = dict(centroid_options or {})
    weight_kwargs = dict(channel_weight_options or {})
    if "fields" in centroid_kwargs or "fields" in weight_kwargs:
        raise ValueError("field state is not accepted by this fixed no-field composition")
    if "norm_type" in weight_kwargs:
        norm_type = weight_kwargs.pop("norm_type")
    else:
        norm_type = "Diagonal Covariance"
    if "outlier_distance" not in weight_kwargs:
        raise ValueError("channel_weight_options must include outlier_distance")

    memberships = composition.loop_entry.initial_membership.class_memberships
    segmentation = composition.thinning.run.segmentation
    initial_lesion = np.asarray(lesion_centroids, dtype=np.float32)
    initial_black = np.asarray(black_hole_centroids, dtype=np.float32)
    values = np.asarray(images, dtype=np.float32)
    types = tuple(object_types)

    updated = compute_lesion_toads_centroids(
        values,
        memberships,
        segmentation,
        class_centroids,
        types,
        lesion_membership=lesion_membership,
        black_hole_membership=black_hole_membership,
        lesion_centroids=initial_lesion,
        black_hole_centroids=initial_black,
        **centroid_kwargs,
    )
    if updated.lesion_centroids is None:
        raise ValueError("post-thinning refresh requires an ordinary or optimized lesion class")

    # The selected wrapper has linkedBlackHoleCentroid=false and computes the
    # direct lesion weighted mean on this active branch. Match that branch with
    # the small source-specific correction used by the standalone refresh.
    direct_lesion = _source_direct_lesion_centroid(
        values,
        memberships,
        segmentation,
        lesion_membership,
        black_hole_membership,
        initial_lesion,
        types,
        fields=None,
        lesion_power_exponents=centroid_kwargs.get("lesion_power_exponents"),
        optimized_auxiliary_power_exponents=centroid_kwargs.get(
            "optimized_auxiliary_power_exponents"
        ),
        use_lesion_weight=use_lesion_weight,
        flair_channel=flair_channel,
    )
    updated = LesionToadsAuxiliaryCentroids(
        updated.class_centroids, direct_lesion, updated.black_hole_centroids
    )
    updated_black = (
        initial_black if updated.black_hole_centroids is None
        else updated.black_hole_centroids
    )
    variance = compute_lesion_toads_variances(
        values,
        memberships,
        segmentation,
        updated.class_centroids,
        lesion_membership,
        black_hole_membership,
        updated.lesion_centroids,
        updated_black,
        types,
    )
    weights = compute_lesion_toads_channel_weights(
        values,
        memberships,
        updated.class_centroids,
        lesion_membership,
        black_hole_membership,
        outlier_membership,
        updated.lesion_centroids,
        updated_black,
        types,
        norm_type=norm_type,
        **weight_kwargs,
    )
    return LesionPostThinningStatistics(
        updated,
        variance,
        weights,
        composition.stage_order + (
            "computeCentroids", "computeVariances", "computeChannelWeights"
        ),
    )
