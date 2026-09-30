"""Compose the next LesionTOADS iteration from a continuing branch result.

Paper context: Shiee et al. (2010), ``docs/papers/lesion-toads.md`` describes
iterative lesion-aware statistical/topological segmentation. This composes a
bounded selected-source transition with fixed inputs, not convergence or a
complete LesionTOADS workflow.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_growth import LesionToadsThinningPass, propagate_lesion_toads_thinning
from .lesion_post_registration import LesionToadsPostRegistrationBranch
from .lesion_relation_aware_membership import (
    LesionRelationAwareMembership,
    run_lesion_toads_relation_aware_membership_pass,
)
from .lesion_statistics_membership_bridge import (
    LesionToadsStatisticsMembershipBridge,
    run_lesion_toads_statistics_membership_bridge,
)


@dataclass(frozen=True)
class LesionToadsContinuingIteration:
    """One continuing iteration through the next energy-ratio calculation."""

    thinning: LesionToadsThinningPass
    statistics_energy: LesionToadsStatisticsMembershipBridge
    relation_membership: LesionRelationAwareMembership
    energy_ratio: np.float32
    stage_order: tuple[str, ...]


def run_lesion_toads_continuing_iteration(
    previous_branch: LesionToadsPostRegistrationBranch,
    images: ArrayLike,
    transformed_shape: ArrayLike,
    class_shape_values: ArrayLike,
    *,
    volume_ratio: float,
    thinning_limit: float,
    bridge_options: Mapping[str, Any],
    use_lesion_max: bool = False,
) -> LesionToadsContinuingIteration:
    """Continue a standard-mode iteration through thinning and memberships.

    The prior branch must have continued and refreshed its relation list. Its
    growth outputs become the next thinning state. The selected wrapper does
    not rebuild that relation list after thinning; the same relation map is
    therefore passed into the following membership update, where direct vs.
    prior is selected per voxel/class. Correction fields, atlas shapes, and
    other spatial inputs in ``bridge_options`` are fixed for this slice.
    """
    if not isinstance(previous_branch, LesionToadsPostRegistrationBranch):
        raise TypeError("previous_branch must be a LesionToadsPostRegistrationBranch")
    if previous_branch.stop or previous_branch.relations is None:
        raise ValueError("the previous wrapper branch must be continuing with refreshed relations")
    if previous_branch.growing is None:
        raise ValueError("this composition requires the standard-mode continuing growth state")
    if not isinstance(bridge_options, Mapping):
        raise TypeError("bridge_options must be a mapping of fixed statistics/spatial inputs")
    prior_result = previous_branch.membership_result
    prior_membership = prior_result.membership
    types = tuple(str(value).lower() for value in bridge_options["object_types"])
    if types != ("object", "lesion", "blackhole", "outmask"):
        raise ValueError("the composed bridge currently supports OBJECT/LESION/BLACKHOLE/OUTMASK")
    prior_stats = prior_result.statistics

    thinning = propagate_lesion_toads_thinning(
        previous_branch.segmentation,
        prior_membership.class_memberships,
        types,
        prior_membership.lesion_membership,
        prior_membership.black_hole_membership,
        transformed_shape,
        previous_branch.growing.skeleton,
        previous_branch.growing.available,
        previous_branch.growing.ordering,
        lesion_class=types.index("lesion"),
        volume_ratio=volume_ratio,
        limit=thinning_limit,
        use_lesion_max=use_lesion_max,
    )

    options = dict(bridge_options)
    options.pop("object_types")
    # The standard bridge supplies the already source-matched statistics and
    # raw field-aware energy construction. Its static-relation membership map
    # is deliberately replaced below by per-voxel relation dispatch.
    options.pop("relation_types", None)
    statistics_energy = run_lesion_toads_statistics_membership_bridge(
        images,
        prior_membership.class_memberships,
        thinning.run.segmentation,
        prior_stats.centroids.class_centroids,
        prior_membership.lesion_membership,
        prior_membership.black_hole_membership,
        prior_membership.outlier_membership,
        prior_stats.centroids.lesion_centroids,
        prior_stats.centroids.black_hole_centroids,
        types,
        relation_types=("direct",) * len(types),
        previous_energy=float(prior_membership.objective),
        **options,
    )
    base_membership = statistics_energy.membership
    class_data = np.zeros_like(base_membership.class_memberships, dtype=np.float32)
    class_data[..., :3] = statistics_energy.data_energies.class_energy
    class_data[..., 3] = statistics_energy.outmask_class_energy
    dynamic = run_lesion_toads_relation_aware_membership_pass(
        prior_membership.class_memberships,
        prior_membership.lesion_membership,
        prior_membership.black_hole_membership,
        prior_membership.outlier_membership,
        class_data,
        statistics_energy.data_energies.lesion_energy,
        statistics_energy.data_energies.black_hole_energy,
        statistics_energy.outmask_data_energy,
        class_shape_values,
        np.asarray(options["class_smoothing_priors"], dtype=np.float32),
        tuple(value.upper() for value in types),
        previous_branch.relations,
        smoothing=float(options["smoothing"]),
        shape_priors=options["shape_priors"],
        class_similarities=statistics_energy.similarity.class_similarity,
        lesion_similarities=statistics_energy.similarity.lesion_similarity,
        black_hole_similarities=statistics_energy.similarity.black_hole_similarity,
        prior_coefficient=float(options["prior_coefficient"]),
    )
    old_energy = np.float32(prior_membership.objective)
    next_energy = np.float32(dynamic.membership.objective)
    numerator = np.float32(np.float32(2.0) * np.float32(old_energy - next_energy))
    denominator = np.float32(old_energy + next_energy)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.float32(np.abs(np.float32(numerator / denominator)))
    completed_statistics_energy = replace(
        statistics_energy,
        membership=dynamic.membership,
        energy_ratio=ratio,
    )
    return LesionToadsContinuingIteration(
        thinning=thinning,
        statistics_energy=completed_statistics_energy,
        relation_membership=dynamic,
        energy_ratio=ratio,
        stage_order=(
            "propagateThinning",
            "computeCentroids",
            "computeVariances",
            "computeChannelWeights",
            "computeMembershipsDataEnergy",
            "computeMembershipsRelationAwareSweep",
            "energyRatio",
        ),
    )
