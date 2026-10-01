"""Refresh the selected class statistics after a LesionTOADS thinning pass.

This seam reuses the general mixed-category source ports for centroids,
variances, and channel weights. It deliberately does not route through the
four-class OBJECT/LESION/BLACKHOLE/OUTMASK membership bridge: MaCRUISEplus uses
its 13-class profile and must retain each selected class category.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_growth import LesionToadsThinningPass
from .lesion_statistics import (
    LesionToadsAuxiliaryCentroids,
    LesionToadsVariances,
    compute_lesion_toads_centroids,
    compute_lesion_toads_channel_weights,
    compute_lesion_toads_variances,
)
from .lesion_statistics_refresh import _source_direct_lesion_centroid


@dataclass(frozen=True)
class LesionToadsNextStatisticsRefresh:
    centroids: LesionToadsAuxiliaryCentroids
    variances: LesionToadsVariances
    channel_weights: NDArray[np.float32]


@dataclass(frozen=True)
class LesionNextIterationStatistics:
    """Statistics refreshed from the segmentation after next-loop thinning."""

    statistics: LesionToadsNextStatisticsRefresh
    stage_order: tuple[str, ...]


def refresh_statistics_after_next_iteration_thinning(
    thinning: LesionToadsThinningPass,
    images: ArrayLike,
    memberships: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    outlier_membership: ArrayLike,
    initial_class_centroids: ArrayLike,
    initial_lesion_centroids: ArrayLike,
    initial_black_hole_centroids: ArrayLike,
    object_types: Sequence[str],
    *,
    fields: ArrayLike | None = None,
    centroid_options: Mapping[str, Any] | None = None,
    channel_weight_options: Mapping[str, Any],
    use_lesion_weight: bool = False,
    flair_channel: int | None = None,
) -> LesionNextIterationStatistics:
    """Run ``computeCentroids`` → ``computeVariances`` → ``computeChannelWeights``.

    The freshly thinned segmentation supplies class labels. Image channels,
    memberships, auxiliary maps, initialized statistics, and any corrected
    fields are explicit fixed inputs. The wrapper's class types are passed
    unchanged to the source-matched generic statistics kernels, so categories
    such as MASK and optimized tissue classes are not collapsed into the
    four-class bridge profile. This stops before the following membership
    update, whose 13-class source composition remains a separate task.
    """
    if not isinstance(thinning, LesionToadsThinningPass):
        raise TypeError("thinning must be a LesionToadsThinningPass")
    values = np.asarray(images, dtype=np.float32)
    mems = np.asarray(memberships, dtype=np.float32)
    labels = np.asarray(thinning.run.segmentation)
    types = tuple(str(value).lower() for value in object_types)
    if values.ndim != 4 or labels.shape != values.shape[:3]:
        raise ValueError("thinned segmentation must match the 3-D image grid")
    if mems.shape[:3] != labels.shape or mems.ndim != 4:
        raise ValueError("class memberships must match the thinned image grid")
    if len(types) != mems.shape[-1]:
        raise ValueError("object_types must match the membership class count")
    centroid_options = dict(centroid_options or {})
    weight_options = dict(channel_weight_options)
    if "fields" in centroid_options or "fields" in weight_options:
        raise ValueError("pass corrected fields only through the fields argument")
    norm_type = weight_options.pop("norm_type", "Diagonal Covariance")
    if "outlier_distance" not in weight_options:
        raise ValueError("channel_weight_options must include outlier_distance")
    lesion_map = np.asarray(lesion_membership, dtype=np.float32)
    black_map = np.asarray(black_hole_membership, dtype=np.float32)
    outlier_map = np.asarray(outlier_membership, dtype=np.float32)
    initial_lesion = np.asarray(initial_lesion_centroids, dtype=np.float32)
    initial_black = np.asarray(initial_black_hole_centroids, dtype=np.float32)
    updated = compute_lesion_toads_centroids(
        values,
        mems,
        labels,
        initial_class_centroids,
        types,
        fields=fields,
        lesion_membership=lesion_map,
        black_hole_membership=black_map,
        lesion_centroids=initial_lesion,
        black_hole_centroids=initial_black,
        **centroid_options,
    )
    if updated.lesion_centroids is None:
        raise ValueError("next-iteration refresh requires a LESION class")
    flair = flair_channel
    if flair is None:
        modalities = tuple(centroid_options.get("modalities") or ())
        flair = max((i for i, name in enumerate(modalities) if name == "FLAIR"), default=-1)
    updated_lesion = _source_direct_lesion_centroid(
        values,
        mems,
        labels,
        lesion_map,
        black_map,
        initial_lesion,
        types,
        fields=fields,
        lesion_power_exponents=centroid_options.get("lesion_power_exponents"),
        optimized_auxiliary_power_exponents=centroid_options.get(
            "optimized_auxiliary_power_exponents"
        ),
        use_lesion_weight=use_lesion_weight,
        flair_channel=flair,
    )
    updated = LesionToadsAuxiliaryCentroids(
        updated.class_centroids, updated_lesion, updated.black_hole_centroids
    )
    black_centroid = (
        initial_black
        if updated.black_hole_centroids is None
        else updated.black_hole_centroids
    )
    variances = compute_lesion_toads_variances(
        values,
        mems,
        labels,
        updated.class_centroids,
        lesion_map,
        black_map,
        updated.lesion_centroids,
        black_centroid,
        types,
        fields=fields,
    )
    weights = compute_lesion_toads_channel_weights(
        values,
        mems,
        updated.class_centroids,
        lesion_map,
        black_map,
        outlier_map,
        updated.lesion_centroids,
        black_centroid,
        types,
        fields=fields,
        norm_type=norm_type,
        **weight_options,
    )
    stats = LesionToadsNextStatisticsRefresh(
        updated, variances, weights
    )
    return LesionNextIterationStatistics(
        stats,
        (
            "propagateThinning",
            "computeCentroids",
            "computeVariances",
            "computeChannelWeights",
        ),
    )
