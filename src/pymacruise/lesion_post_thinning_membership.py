"""Compose post-thinning statistics into the next standard-13 membership pass.

This fixed-input continuation follows the selected wrapper order after the
initial thinning: refresh centroids, variances, and channel weights, then run
the next ``computeMemberships`` volume pass with ``GMflag=true``. The shape,
relation, and distance maps are supplied by the caller and remain fixed; this
adapter does not synthesize absent classes or generate those maps.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any
from collections.abc import Mapping

from numpy.typing import ArrayLike

from .lesion_initial_thinning_composition import LesionInitialThinningComposition
from .lesion_post_thinning_statistics import (
    LesionPostThinningStatistics,
    refresh_statistics_after_initial_thinning,
)
from .lesion_standard13_volume import Standard13VolumeResult, run_lesion_toads_standard13_volume


@dataclass(frozen=True)
class LesionPostThinningMembershipUpdate:
    statistics: LesionPostThinningStatistics
    membership: Standard13VolumeResult
    stage_order: tuple[str, ...]


def run_post_thinning_standard13_membership_update(
    composition: LesionInitialThinningComposition,
    images: ArrayLike,
    initial_class_centroids: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    outlier_membership: ArrayLike,
    initial_lesion_centroid: ArrayLike,
    initial_black_hole_centroid: ArrayLike,
    shape_values: ArrayLike,
    shape_priors: ArrayLike,
    smoothing_priors: ArrayLike,
    direct_relations: ArrayLike,
    initial_lesion_relation: ArrayLike,
    ventricle_distances: ArrayLike,
    gray_matter_distances: ArrayLike,
    brainstem_distances: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    *,
    centroid_options: Mapping[str, Any] | None = None,
    channel_weight_options: Mapping[str, Any],
    outlier_energy: float,
    lesion_power_exponent: float,
    use_lesion_weight: bool = False,
    flair_channel: int = -1,
) -> LesionPostThinningMembershipUpdate:
    """Refresh source statistics and feed them into the next GM-enabled pass."""
    stats = refresh_statistics_after_initial_thinning(
        composition,
        images,
        initial_class_centroids,
        lesion_membership,
        black_hole_membership,
        outlier_membership,
        initial_lesion_centroid,
        initial_black_hole_centroid,
        object_types,
        centroid_options=centroid_options,
        channel_weight_options=channel_weight_options,
        use_lesion_weight=use_lesion_weight,
        flair_channel=flair_channel,
    )
    lesion_centroids = stats.centroids.lesion_centroids
    if lesion_centroids is None:
        raise ValueError("the standard-13 update requires a refreshed LESION centroid")
    lesion_stddev = stats.variances.lesion_standard_deviation
    if lesion_stddev is None:
        raise ValueError("the standard-13 update requires a refreshed LESION deviation")

    initial = composition.loop_entry.initial_membership
    membership = run_lesion_toads_standard13_volume(
        images=images,
        class_centroids=stats.centroids.class_centroids,
        class_stddev=stats.variances.class_standard_deviation,
        lesion_centroid=lesion_centroids,
        lesion_stddev=lesion_stddev,
        channel_weights=stats.channel_weights,
        shape_values=shape_values,
        shape_priors=shape_priors,
        smoothing_priors=smoothing_priors,
        direct_relations=direct_relations,
        initial_class_memberships=initial.class_memberships,
        initial_lesion_membership=initial.lesion_membership,
        initial_black_hole_membership=initial.black_hole_membership,
        initial_outlier_membership=initial.outlier_membership,
        initial_lesion_relation=initial.lesion_relation,
        ventricle_distances=ventricle_distances,
        gray_matter_distances=gray_matter_distances,
        brainstem_distances=brainstem_distances,
        outlier_energy=outlier_energy,
        lesion_power_exponent=lesion_power_exponent,
        gm_flag=True,
        source_zero_variance_nan=True,
    )
    return LesionPostThinningMembershipUpdate(
        stats,
        membership,
        stats.stage_order + ("computeMemberships:GMflag=true",),
    )
