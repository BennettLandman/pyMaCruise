"""Run the selected 13-class GM-enabled membership and energy update."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_energy_ratio import LesionEnergyRatio, lesion_toads_energy_ratio
from .lesion_next_iteration_statistics import LesionNextIterationStatistics
from .lesion_standard13_volume import Standard13VolumeResult, run_lesion_toads_standard13_volume


@dataclass(frozen=True)
class LesionNextIterationMembership:
    """Updated standard-13 memberships and source-ordered energy ratio."""

    membership: Standard13VolumeResult
    energy_ratio: LesionEnergyRatio
    stage_order: tuple[str, ...]


def run_lesion_toads_next_gm_membership_update(
    statistics: LesionNextIterationStatistics,
    images: ArrayLike,
    shape_values: ArrayLike,
    shape_priors: ArrayLike,
    smoothing_priors: ArrayLike,
    direct_relations: ArrayLike,
    initial_class_memberships: ArrayLike,
    initial_lesion_membership: ArrayLike,
    initial_black_hole_membership: ArrayLike,
    initial_outlier_membership: ArrayLike,
    initial_lesion_relation: ArrayLike,
    ventricle_distances: ArrayLike,
    gray_matter_distances: ArrayLike,
    brainstem_distances: ArrayLike,
    *,
    previous_energy: float,
    outlier_energy: float,
    lesion_power_exponent: float,
    spread: float = 0.5,
    source_zero_variance_nan: bool = False,
) -> LesionNextIterationMembership:
    """Call the fixed MaCRUISEplus 13-class volume kernel with ``GMflag=true``.

    The direct/prior class flags, transformed shape samples, lesion relation,
    auxiliary maps, and distance maps belong to the prior wrapper state and
    are explicit fixed inputs for this membership sweep. This calls the
    standard-13 source-order kernel rather than the separate four-class
    relation-aware bridge. The returned energy is compared with the supplied
    previous objective using the selected float32 wrapper expression.
    """
    if not isinstance(statistics, LesionNextIterationStatistics):
        raise TypeError("statistics must be the preceding 13-class statistics refresh")
    refreshed = statistics.statistics
    membership = run_lesion_toads_standard13_volume(
        images,
        refreshed.centroids.class_centroids,
        refreshed.variances.class_standard_deviation,
        refreshed.centroids.lesion_centroids,
        refreshed.variances.lesion_standard_deviation,
        refreshed.channel_weights,
        shape_values,
        shape_priors,
        smoothing_priors,
        direct_relations,
        initial_class_memberships,
        initial_lesion_membership,
        initial_black_hole_membership,
        initial_outlier_membership,
        initial_lesion_relation,
        ventricle_distances,
        gray_matter_distances,
        brainstem_distances,
        outlier_energy=outlier_energy,
        lesion_power_exponent=lesion_power_exponent,
        gm_flag=True,
        spread=spread,
        source_zero_variance_nan=source_zero_variance_nan,
    )
    ratio = lesion_toads_energy_ratio(previous_energy, membership.objective)
    return LesionNextIterationMembership(
        membership,
        ratio,
        statistics.stage_order + ("computeMemberships:GMflag=true", "energyRatio"),
    )
