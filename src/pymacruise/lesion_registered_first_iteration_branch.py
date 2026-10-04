"""Carry a registered LesionTOADS first loop update through branch control."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from numpy.typing import ArrayLike

from .lesion_next_iteration_branch import (
    LesionNextIterationBranch,
    continue_after_next_iteration_registration,
)
from .lesion_next_iteration_registration import (
    LesionNextIterationRegistration,
    register_after_next_gm_membership_update,
)
from .lesion_registered_first_iteration import LesionRegisteredFirstIteration
from .lesion_registered_initial_membership import LesionRegisteredInitialMembership


@dataclass(frozen=True)
class LesionRegisteredFirstIterationBranch:
    """First loop update carried through registration and wrapper control."""

    transition: LesionRegisteredFirstIteration
    registration: LesionNextIterationRegistration
    branch: LesionNextIterationBranch
    stage_order: tuple[str, ...]


def run_lesion_toads_registered_first_iteration_branch(
    initial: LesionRegisteredInitialMembership,
    transition: LesionRegisteredFirstIteration,
    atlas_shape_priors: Sequence[ArrayLike],
    initial_transform: ArrayLike,
    *,
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
    align_shapes: bool = True,
    init_align_iter: int = 10,
    has_lesion: bool = True,
    use_lesion_max: bool = False,
) -> LesionRegisteredFirstIterationBranch:
    """Run first-loop membership → atlas registration → source branch.

    The supplied ``initial`` state contributes the first thinning result and
    constructor classification/skeleton/frontier arrays. The transition
    contributes the post-thinning membership and energy ratio. Registration
    is executed before the selected convergence and iteration checks. Any
    continuing growth/relation refresh or terminal classification follows
    the selected branch. This adapter ends after that one branch and does not
    initialize another iteration.
    """
    if not isinstance(initial, LesionRegisteredInitialMembership):
        raise TypeError("initial must be a registered initial-membership result")
    if not isinstance(transition, LesionRegisteredFirstIteration):
        raise TypeError("transition must be a registered first-iteration result")
    registration = register_after_next_gm_membership_update(
        transition.membership,
        atlas_shape_priors,
        initial_transform,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        align_shapes=align_shapes,
        init_align_iter=init_align_iter,
    )
    thinning = initial.first_thinning.thinning.run
    branch = continue_after_next_iteration_registration(
        registration,
        segmentation=thinning.segmentation,
        skeleton=thinning.skeleton,
        ordering=thinning.ordering,
        class_volume=thinning.class_volume,
        object_types=object_types,
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
    stages = (
        *initial.stage_order,
        *transition.stage_order[1:],
        *branch.stage_order[len(transition.stage_order):],
    )
    return LesionRegisteredFirstIterationBranch(
        transition,
        registration,
        branch,
        stages,
    )
