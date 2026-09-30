"""Compose the LesionTOADS wrapper branch after a continuing iteration.

This adapter consumes the real P39zzzzv result, registers its updated
memberships, then applies the selected convergence/iteration-limit branch.
It returns the existing post-registration branch state so a continuing result
can be passed directly to :func:`run_lesion_toads_continuing_iteration`.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike

from .lesion_continuing_iteration import LesionToadsContinuingIteration
from .lesion_post_registration import (
    LesionToadsPostRegistrationBranch,
    lesion_toads_post_registration_branch,
)
from .lesion_registration_handoff import register_lesion_membership_handoff


def run_lesion_toads_continuing_post_membership_branch(
    previous_branch: LesionToadsPostRegistrationBranch,
    continuing_iteration: LesionToadsContinuingIteration,
    shape_priors: Sequence[ArrayLike],
    *,
    previous_energy: float,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    object_types: Sequence[str],
    classification: ArrayLike,
    template_labels: ArrayLike,
    lesion_class: int,
    iteration: int,
    n_iter_max: int,
    max_distance: float,
    algorithm_mode: str,
    relation_limit: int,
    has_lesion: bool = True,
    use_lesion_max: bool = False,
    align_shapes: bool = True,
    init_align_iter: int = 10,
) -> LesionToadsPostRegistrationBranch:
    """Run registration and the wrapper branch after P39zzzzv.

    ``previous_energy`` and the iteration controls are explicit because the
    selected wrapper carries energy/control state across loop boundaries.
    The supplied previous energy must equal the prior branch's objective;
    this prevents accidentally computing the ratio against P39v's new energy.
    This composition is restricted to P39v's standard-mode OBJECT/LESION/
    BLACKHOLE/OUTMASK path.
    """
    if not isinstance(previous_branch, LesionToadsPostRegistrationBranch):
        raise TypeError("previous_branch must be a LesionToadsPostRegistrationBranch")
    if previous_branch.stop or previous_branch.relations is None or previous_branch.growing is None:
        raise ValueError("previous_branch must be a continuing result with refreshed relations")
    if not isinstance(continuing_iteration, LesionToadsContinuingIteration):
        raise TypeError("continuing_iteration must be a P39zzzzv result")
    if not isinstance(algorithm_mode, str) or algorithm_mode != "standard":
        raise ValueError("P39zzzzv post-membership composition currently supports standard mode only")
    expected_types = ("object", "lesion", "blackhole", "outmask")
    types = tuple(str(value).lower() for value in object_types)
    if types != expected_types:
        raise ValueError("object_types must be OBJECT/LESION/BLACKHOLE/OUTMASK in source order")
    previous = np.float32(previous_energy)
    source_previous = np.float32(previous_branch.membership_result.membership.objective)
    if not np.isfinite(previous) or previous.view(np.uint32) != source_previous.view(np.uint32):
        raise ValueError("previous_energy must match the prior branch's float32 objective")
    if continuing_iteration.statistics_energy.membership is not continuing_iteration.relation_membership.membership:
        raise ValueError("P39v statistics and relation-aware outputs do not share membership state")

    registration = register_lesion_membership_handoff(
        continuing_iteration.statistics_energy,
        shape_priors,
        previous_branch.registration_handoff.transform,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        align_shapes=align_shapes,
        init_align_iter=init_align_iter,
    )
    thinning_run = continuing_iteration.thinning.run
    branch = lesion_toads_post_registration_branch(
        continuing_iteration.statistics_energy,
        registration,
        previous_energy=float(previous),
        segmentation=thinning_run.segmentation,
        skeleton=thinning_run.skeleton,
        ordering=thinning_run.ordering,
        class_volume=thinning_run.class_volume,
        object_types=types,
        classification=classification,
        template_labels=template_labels,
        lesion_class=lesion_class,
        iteration=iteration,
        n_iter_max=n_iter_max,
        max_distance=max_distance,
        algorithm_mode=algorithm_mode,
        relation_limit=relation_limit,
        has_lesion=has_lesion,
        use_lesion_max=use_lesion_max,
    )
    # P39v has already emitted the ratio before registration. The combined
    # trace follows wrapper order and keeps this branch directly reusable as
    # the next P39v starting state.
    stages = list(continuing_iteration.stage_order)
    stages.extend(("registerShapes", "convergenceCheck", "iterationLimitCheck"))
    if branch.run_growing:
        stages.append("propagateGrowing")
    if branch.stop:
        stages.extend(("updateClassificationFromLabels", "updateLesionClassifcation"))
    else:
        stages.append("computeRelations")
    return replace(branch, stage_order=tuple(stages))
