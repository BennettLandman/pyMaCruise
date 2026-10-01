"""Cross-check membership and competition adapters on folded synthetic state.

This diagnostic composition combines outputs from separate selected modes: the
default ``linear_fast_marching`` thinning/statistics/membership path and the
alternative ``narrow_band_competition`` queue path. The selected wrapper does
not execute competition after this membership/ratio update, so this helper is
not an end-to-end or source-order workflow. It checks fixed-input adapter
interoperability only; initialization, registration, convergence, and the
complete LesionTOADS workflow remain outside scope.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from numpy.typing import ArrayLike

from .lesion_competition import LesionToadsCompetitionPass, propagate_lesion_toads_competition
from .lesion_energy_ratio import LesionEnergyRatio, lesion_toads_energy_ratio
from .lesion_initial_thinning_composition import LesionInitialThinningComposition
from .lesion_post_thinning_membership import (
    LesionPostThinningMembershipUpdate,
    run_post_thinning_standard13_membership_update,
)


@dataclass(frozen=True)
class LesionPostThinningCompetition:
    """Selected continuation and all exposed stage results."""

    membership_update: LesionPostThinningMembershipUpdate
    energy: LesionEnergyRatio
    competition: LesionToadsCompetitionPass
    stage_order: tuple[str, ...]


def run_post_thinning_competition_pass(
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
    lesion_class: int,
    use_lesion_weight: bool = False,
    flair_channel: int = -1,
) -> LesionPostThinningCompetition:
    """Run selected post-thinning refresh, GM pass, energy ratio, and competition.

    Competition consumes the GMflag=true class/auxiliary memberships and the
    segmentation, ordering, and class-volume states emitted by first thinning.
    The prior objective is the GMflag=false initial membership objective.
    Statistical NaNs produced by the selected zero-standard-deviation source
    arithmetic are passed through without replacing or masking them.
    """
    if not isinstance(composition, LesionInitialThinningComposition):
        raise TypeError("composition must be a LesionInitialThinningComposition")
    update = run_post_thinning_standard13_membership_update(
        composition,
        images,
        initial_class_centroids,
        lesion_membership,
        black_hole_membership,
        outlier_membership,
        initial_lesion_centroid,
        initial_black_hole_centroid,
        shape_values,
        shape_priors,
        smoothing_priors,
        direct_relations,
        initial_lesion_relation,
        ventricle_distances,
        gray_matter_distances,
        brainstem_distances,
        object_types,
        centroid_options=centroid_options,
        channel_weight_options=channel_weight_options,
        outlier_energy=outlier_energy,
        lesion_power_exponent=lesion_power_exponent,
        use_lesion_weight=use_lesion_weight,
        flair_channel=flair_channel,
    )
    initial = composition.loop_entry.initial_membership
    energy = lesion_toads_energy_ratio(
        float(initial.objective), float(update.membership.objective)
    )
    thinning = composition.thinning.run
    competition = propagate_lesion_toads_competition(
        thinning.segmentation,
        update.membership.class_memberships,
        update.membership.lesion_membership,
        update.membership.black_hole_membership,
        thinning.ordering,
        thinning.class_volume,
        lesion_class=lesion_class,
        object_types=tuple(object_types),
    )
    return LesionPostThinningCompetition(
        update,
        energy,
        competition,
        update.stage_order + ("energyRatio", "initCompetitionLabels", "propagateCompetition"),
    )
