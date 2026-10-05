"""Continue registered LesionTOADS through its first loop membership update.

This composes the first thinning handoff from the registered initial
membership stage with the selected post-thinning statistics, GM-enabled
membership, and energy-ratio operations. It does not perform registration,
field correction, convergence branching, growth, or a complete run.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike

from .lesion_next_iteration_membership import (
    LesionNextIterationMembership,
    run_lesion_toads_next_gm_membership_update,
)
from .lesion_next_iteration_statistics import (
    LesionNextIterationStatistics,
    refresh_statistics_after_next_iteration_thinning,
)
from .lesion_registered_initial_membership import (
    LesionInitialDistanceMaps,
    LesionRegisteredInitialMembership,
)
from .lesion_distance_refresh import refresh_lesion_toads_distance_maps
from .lesion_standard13_voxel import STANDARD13_TYPES


@dataclass(frozen=True)
class LesionRegisteredFirstIteration:
    """Post-thinning statistics, next memberships, and energy ratio."""

    statistics: LesionNextIterationStatistics
    membership: LesionNextIterationMembership
    stage_order: tuple[str, ...]
    refreshed_distance_maps: LesionInitialDistanceMaps | None = None


def run_lesion_toads_registered_first_iteration(
    initial: LesionRegisteredInitialMembership,
    images: ArrayLike,
    initial_class_centroids: ArrayLike,
    initial_lesion_centroids: ArrayLike,
    initial_black_hole_centroids: ArrayLike,
    transformed_shape_values: ArrayLike,
    shape_priors: ArrayLike,
    smoothing_priors: ArrayLike,
    *,
    centroid_options: Mapping[str, Any],
    channel_weight_options: Mapping[str, Any],
    outlier_energy: float,
    lesion_power_exponent: float,
    spread: float = 1.0,
    source_zero_variance_nan: bool | None = None,
    atlas_names: Sequence[str] | None = None,
    maximum_ventricle_distance: int = 2,
    maximum_gm_distance: int = 3,
    maximum_brainstem_distance: int = 25,
    has_lesion: bool = True,
) -> LesionRegisteredFirstIteration:
    """Run initial thinning → refreshed statistics → next memberships/ratio.

    Spatial atlas relations and membership maps are carried from ``initial``.
    When atlas names are supplied, source GM-enabled distance maps are rebuilt
    from the post-thinning segmentation before membership; omitting names
    retains the fixed constructor maps for reference callers. Class and
    auxiliary centroids are the values carried into the registered recrop.
    Corrected fields are deliberately not synthesized: the selected
    MaCRUISEplus profile disables that optional correction.
    """
    if not isinstance(initial, LesionRegisteredInitialMembership):
        raise TypeError("initial must be a registered initial-membership result")
    if not isinstance(centroid_options, Mapping):
        raise TypeError("centroid_options must be a mapping")
    if not isinstance(channel_weight_options, Mapping):
        raise TypeError("channel_weight_options must be a mapping")

    previous = initial.membership
    statistics = refresh_statistics_after_next_iteration_thinning(
        initial.first_thinning.thinning,
        images,
        previous.class_memberships,
        previous.lesion_membership,
        previous.black_hole_membership,
        previous.outlier_membership,
        initial_class_centroids,
        initial_lesion_centroids,
        initial_black_hole_centroids,
        STANDARD13_TYPES,
        centroid_options=centroid_options,
        channel_weight_options=channel_weight_options,
        use_lesion_weight=True,
        flair_channel=1,
    )
    direct_relations = initial.relations.relation_membership[
        initial.relations.relation_map
    ]
    class_deviations = statistics.statistics.variances.class_standard_deviation
    lesion_deviations = statistics.statistics.variances.lesion_standard_deviation
    detected_zero_variance_nan = bool(
        np.any(class_deviations == np.float32(0.0))
        or np.any(lesion_deviations == np.float32(0.0))
    )
    if source_zero_variance_nan is None:
        source_zero_variance_nan = detected_zero_variance_nan
    elif not isinstance(source_zero_variance_nan, (bool, np.bool_)):
        raise TypeError("source_zero_variance_nan must be boolean or None")
    if not isinstance(has_lesion, (bool, np.bool_)):
        raise TypeError("has_lesion must be boolean")
    refreshed = (
        refresh_lesion_toads_distance_maps(
            initial.first_thinning.thinning.run.segmentation,
            atlas_names,
            maximum_ventricle_distance=maximum_ventricle_distance,
            maximum_gm_distance=maximum_gm_distance,
            maximum_brainstem_distance=maximum_brainstem_distance,
            previous_maps=initial.distance_maps,
        )
        if atlas_names is not None and has_lesion else None
    )
    distance_maps = initial.distance_maps if refreshed is None else refreshed.maps
    next_membership = run_lesion_toads_next_gm_membership_update(
        statistics,
        images,
        transformed_shape_values,
        shape_priors,
        smoothing_priors,
        direct_relations,
        previous.class_memberships,
        previous.lesion_membership,
        previous.black_hole_membership,
        previous.outlier_membership,
        previous.lesion_relation,
        distance_maps.ventricle,
        distance_maps.gray_matter,
        distance_maps.brainstem,
        previous_energy=float(previous.objective),
        outlier_energy=outlier_energy,
        lesion_power_exponent=lesion_power_exponent,
        spread=spread,
        source_zero_variance_nan=source_zero_variance_nan,
    )
    return LesionRegisteredFirstIteration(
        statistics,
        next_membership,
        (
            *statistics.stage_order,
            *((refreshed.stage_order) if refreshed else ()),
            *next_membership.stage_order[len(statistics.stage_order):],
        ),
        None if refreshed is None else refreshed.maps,
    )
