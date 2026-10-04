"""Refresh selected LesionTOADS statistics and update a descriptor membership."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from numpy.typing import ArrayLike

from .lesion_descriptor_iteration_slice import (
    LesionToadsDescriptorIterationSlice,
    run_lesion_toads_descriptor_iteration_slice,
)
from .lesion_descriptor_membership import map_lesion_descriptor_membership_types
from .lesion_similarity import LesionToadsSimilarity, compute_lesion_toads_similarity
from .lesion_statistics_refresh import LesionToadsStatisticsRefresh, refresh_lesion_toads_statistics


@dataclass(frozen=True)
class LesionToadsDescriptorStatisticsMembership:
    """Statistics, similarities, membership transition, and source energy ratio."""

    statistics: LesionToadsStatisticsRefresh
    similarity: LesionToadsSimilarity
    transition: LesionToadsDescriptorIterationSlice
    energy_ratio: np.float32
    stage_order: tuple[str, ...]


def run_lesion_toads_descriptor_statistics_membership(
    images: ArrayLike,
    initial_class_memberships: ArrayLike,
    segmentation: ArrayLike,
    initial_class_centroids: ArrayLike,
    initial_lesion_membership: ArrayLike,
    initial_black_hole_membership: ArrayLike,
    initial_outlier_membership: ArrayLike,
    initial_lesion_centroid: ArrayLike,
    initial_black_hole_centroid: ArrayLike,
    descriptor_types: Sequence[str],
    structure_names: Sequence[str],
    *,
    fields: ArrayLike | None,
    centroid_options: Mapping[str, Any],
    channel_weight_options: Mapping[str, Any],
    prior_scale: float,
    initial_black_hole_similarity: ArrayLike,
    smoothing: float,
    smoothing_priors: ArrayLike,
    shape_priors: ArrayLike,
    shape_values: ArrayLike,
    relation_types: Sequence[str],
    prior_coefficient: float,
    outlier_distance: float,
    previous_energy: float,
    use_lesion_weight: bool = False,
    flair_channel: int | None = None,
    t1_channel: int = -1,
    zero: float = 1.0e-10,
    infinity: float = 1.0e10,
) -> LesionToadsDescriptorStatisticsMembership:
    """Refresh fixed descriptor statistics, then compose one membership step.

    This ports the wrapper order after thinning for the selected bundled
    descriptor. Inputs include the already updated class/auxiliary maps and
    hard segmentation. Atlas resampling, topology propagation, correction
    field estimation, distance modulation, and loop stopping are outside this
    function. The selected descriptor does not contain BLACKHOLE, but the
    source keeps an independent black-hole map/centroid whose empty-map update
    is represented by the statistics adapter.
    """
    if not isinstance(use_lesion_weight, (bool, np.bool_)):
        raise ValueError("use_lesion_weight must be boolean")
    if use_lesion_weight:
        raise ValueError("the active descriptor energy adapter does not yet expose FLAIR lesion weighting")
    types = tuple(
        "object" if str(value).lower() == "obj" else str(value).lower()
        for value in descriptor_types
    )
    membership_types = map_lesion_descriptor_membership_types(types)
    if len(structure_names) != len(types) or any(not isinstance(value, str) for value in structure_names):
        raise ValueError("structure_names must provide one name per descriptor class")
    stats = refresh_lesion_toads_statistics(
        images,
        initial_class_memberships,
        segmentation,
        initial_class_centroids,
        initial_lesion_membership,
        initial_black_hole_membership,
        initial_outlier_membership,
        initial_lesion_centroid,
        initial_black_hole_centroid,
        types,
        fields=fields,
        centroid_options=centroid_options,
        channel_weight_options=channel_weight_options,
        use_lesion_weight=False,
        flair_channel=flair_channel,
    )
    class_std = stats.variances.class_standard_deviation
    lesion_std = stats.variances.lesion_standard_deviation
    lesion_center = stats.centroids.lesion_centroids
    black_center = stats.centroids.black_hole_centroids
    if lesion_center is None or black_center is None or lesion_std is None:
        raise ValueError("descriptor statistics require lesion-related auxiliary centroids and deviation")
    similarity = compute_lesion_toads_similarity(
        stats.centroids.class_centroids,
        class_std,
        lesion_center,
        lesion_std,
        stats.channel_weights,
        types,
        structure_names,
        prior_scale=prior_scale,
        t1_channel=t1_channel,
        black_hole_similarity=initial_black_hole_similarity,
    )
    transition = run_lesion_toads_descriptor_iteration_slice(
        images,
        initial_class_memberships,
        types,
        stats.centroids.class_centroids.T,
        lesion_center,
        black_center,
        stats.channel_weights,
        shape_values,
        smoothing_priors,
        shape_priors,
        similarity.class_similarity,
        similarity.lesion_similarity,
        similarity.black_hole_similarity,
        relation_types,
        outlier_distance=outlier_distance,
        smoothing=smoothing,
        prior_coefficient=prior_coefficient,
        initial_lesion_membership=initial_lesion_membership,
        initial_black_hole_membership=initial_black_hole_membership,
        initial_outlier_membership=initial_outlier_membership,
        use_field=fields is not None,
        intensity_fields=fields,
        zero=zero,
        infinity=infinity,
    )
    old = np.float32(previous_energy)
    new = np.float32(transition.membership.objective)
    numerator = np.float32(np.float32(2.0) * np.float32(old - new))
    denominator = np.float32(old + new)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.float32(np.abs(np.float32(numerator / denominator)))
    return LesionToadsDescriptorStatisticsMembership(
        stats,
        similarity,
        transition,
        ratio,
        (
            "computeCentroids",
            "computeVariances",
            "computeChannelWeights",
            "computeCentroidSimilarities",
            "computeActiveClassDataEnergies",
            "computeMemberships",
            "computeEnergyRatio",
        ),
    )
