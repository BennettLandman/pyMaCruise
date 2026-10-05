"""Compose a continuing selected 13-class LesionTOADS update.

The transition joins the source-ordered next thinning, statistics refresh,
GM-enabled distance refresh (when atlas names are supplied), volume membership
update, and energy ratio. It reuses selected VandyCRUISE-backed adapters;
registration, field fitting, post-registration growth, and workflow I/O remain
separate stages.
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
from .lesion_next_iteration_thinning import run_lesion_toads_next_iteration_thinning
from .lesion_post_registration import LesionToadsPostRegistrationBranch
from .lesion_growth import LesionToadsThinningPass
from .lesion_registered_initial_membership import LesionInitialDistanceMaps
from .lesion_distance_refresh import refresh_lesion_toads_distance_maps
from .lesion_inhomogeneity_transition import (
    LesionToadsInhomogeneityTransition,
    update_lesion_toads_exact_inhomogeneity_fields,
)


@dataclass(frozen=True)
class LesionToadsDescriptorContinuingIteration:
    """Outputs from one selected 13-class continuing iteration transition."""

    thinning: LesionToadsThinningPass
    statistics: LesionNextIterationStatistics
    membership: LesionNextIterationMembership
    stage_order: tuple[str, ...]
    refreshed_distance_maps: LesionInitialDistanceMaps | None = None
    inhomogeneity: LesionToadsInhomogeneityTransition | None = None


def run_lesion_toads_descriptor_continuing_iteration(
    previous_branch: LesionToadsPostRegistrationBranch,
    images: ArrayLike,
    transformed_shape_values: ArrayLike,
    *,
    object_types: Sequence[str],
    lesion_class: int,
    volume_ratio: float,
    thinning_limit: float,
    initial_class_centroids: ArrayLike,
    initial_lesion_centroids: ArrayLike,
    initial_black_hole_centroids: ArrayLike,
    shape_priors: ArrayLike,
    smoothing_priors: ArrayLike,
    direct_relations: ArrayLike,
    ventricle_distances: ArrayLike,
    gray_matter_distances: ArrayLike,
    brainstem_distances: ArrayLike,
    centroid_options: Mapping[str, Any] | None,
    channel_weight_options: Mapping[str, Any],
    outlier_energy: float,
    lesion_power_exponent: float,
    fields: ArrayLike | None = None,
    flair_channel: int | None = None,
    spread: float = 0.5,
    source_zero_variance_nan: bool = False,
    use_lesion_max: bool = False,
    atlas_names: Sequence[str] | None = None,
    maximum_ventricle_distance: int = 2,
    maximum_gm_distance: int = 3,
    maximum_brainstem_distance: int = 25,
    has_lesion: bool = True,
    previous_distance_maps: LesionInitialDistanceMaps | None = None,
    correction_polynomial_degree: int | None = None,
    correction_degree: int = 0,
    correction_max_distance: float = 0.0,
) -> LesionToadsDescriptorContinuingIteration:
    """Run thinning → 13-class statistics → GM membership and energy ratio.

    Geometry-dependent inputs (transformed shapes, relation choices and
    distance maps) are the already-updated state handed into this transition.
    The prior class and auxiliary membership maps and objective are read from
    ``previous_branch``. Correction fields and initial statistics are supplied
    explicitly because the branch object does not own those workflow states.
    """
    if not isinstance(previous_branch, LesionToadsPostRegistrationBranch):
        raise TypeError("previous_branch must be a LesionToadsPostRegistrationBranch")
    if previous_branch.stop or previous_branch.relations is None:
        raise ValueError("the previous branch must continue with a refreshed relation state")
    if previous_branch.growing is None:
        raise ValueError("the previous branch must carry its standard-mode growth state")
    if not isinstance(source_zero_variance_nan, (bool, np.bool_)):
        raise TypeError("source_zero_variance_nan must be boolean")
    if not isinstance(has_lesion, (bool, np.bool_)):
        raise TypeError("has_lesion must be boolean")
    values = np.asarray(images, dtype=np.float32)
    membership_update = getattr(previous_branch.membership_result, "membership", None)
    if membership_update is None or not hasattr(membership_update, "class_memberships"):
        raise TypeError("previous branch must carry a 13-class membership result")
    class_memberships = membership_update.class_memberships
    lesion_membership = membership_update.lesion_membership
    black_hole_membership = membership_update.black_hole_membership
    outlier_membership = membership_update.outlier_membership
    lesion_relation = membership_update.lesion_relation
    previous_energy = float(membership_update.objective)

    thinning = run_lesion_toads_next_iteration_thinning(
        previous_branch,
        transformed_shape_values,
        object_types=object_types,
        lesion_class=lesion_class,
        volume_ratio=volume_ratio,
        limit=thinning_limit,
        use_lesion_max=use_lesion_max,
    )
    inhomogeneity = None
    if correction_polynomial_degree is not None:
        if fields is None:
            raise ValueError("correction fields are required when refitting is enabled")
        image_fields = np.asarray(fields, dtype=np.float32)
        if image_fields.shape != values.shape or not np.isfinite(image_fields).all():
            raise ValueError("fields must be finite and match images when refitting is enabled")
        previous_membership = getattr(previous_branch.membership_result, "membership", None)
        if previous_membership is None or not hasattr(previous_membership, "class_memberships"):
            raise TypeError("previous branch lacks memberships required for field refitting")
        inhomogeneity = update_lesion_toads_exact_inhomogeneity_fields(
            values,
            previous_membership.class_memberships,
            thinning.run.segmentation,
            initial_class_centroids,
            tuple(image_fields[..., channel] for channel in range(image_fields.shape[-1])),
            tuple(str(kind).lower() for kind in object_types),
            degree=correction_degree,
            polynomial_degree=correction_polynomial_degree,
            previous_energy_ratio=float(previous_branch.energy_ratio),
            max_distance=correction_max_distance,
        )
        fields = np.stack(inhomogeneity.fields, axis=-1)
    statistics = refresh_statistics_after_next_iteration_thinning(
        thinning,
        values,
        class_memberships,
        lesion_membership,
        black_hole_membership,
        outlier_membership,
        initial_class_centroids,
        initial_lesion_centroids,
        initial_black_hole_centroids,
        object_types,
        fields=fields,
        centroid_options=centroid_options,
        channel_weight_options=channel_weight_options,
        flair_channel=flair_channel,
    )
    refreshed = (
        refresh_lesion_toads_distance_maps(
            thinning.run.segmentation,
            atlas_names,
            maximum_ventricle_distance=maximum_ventricle_distance,
            maximum_gm_distance=maximum_gm_distance,
            maximum_brainstem_distance=maximum_brainstem_distance,
            previous_maps=previous_distance_maps,
        )
        if atlas_names is not None and has_lesion else None
    )
    if refreshed is not None:
        ventricle_distances = refreshed.maps.ventricle
        gray_matter_distances = refreshed.maps.gray_matter
        brainstem_distances = refreshed.maps.brainstem
    membership = run_lesion_toads_next_gm_membership_update(
        statistics,
        values,
        transformed_shape_values,
        shape_priors,
        smoothing_priors,
        direct_relations,
        class_memberships,
        lesion_membership,
        black_hole_membership,
        outlier_membership,
        lesion_relation,
        ventricle_distances,
        gray_matter_distances,
        brainstem_distances,
        previous_energy=previous_energy,
        outlier_energy=outlier_energy,
        lesion_power_exponent=lesion_power_exponent,
        spread=spread,
        source_zero_variance_nan=source_zero_variance_nan,
        intensity_fields=fields,
    )
    return LesionToadsDescriptorContinuingIteration(
        thinning,
        statistics,
        membership,
        (
            *(inhomogeneity.stage_order if inhomogeneity is not None else ()),
            *statistics.stage_order,
            *(refreshed.stage_order if refreshed else ()),
            *membership.stage_order[len(statistics.stage_order):],
        ),
        None if refreshed is None else refreshed.maps,
        inhomogeneity,
    )
