"""Carry a selected 13-class transition through registration and loop control."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from numpy.typing import ArrayLike

from .lesion_descriptor_continuing_iteration import (
    LesionToadsDescriptorContinuingIteration,
)
from .lesion_next_iteration_branch import (
    LesionNextIterationBranch,
    continue_after_next_iteration_registration,
)
from .lesion_next_iteration_registration import (
    LesionNextIterationRegistration,
    register_after_next_gm_membership_update,
)


@dataclass(frozen=True)
class LesionToadsDescriptorRegisteredIteration:
    """Thinning/statistics/membership update plus registration and branch."""

    transition: LesionToadsDescriptorContinuingIteration
    registration: LesionNextIterationRegistration
    control: LesionNextIterationBranch
    stage_order: tuple[str, ...]


def register_and_branch_lesion_toads_descriptor_iteration(
    transition: LesionToadsDescriptorContinuingIteration,
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
    retain_growth_event_states: bool = True,
) -> LesionToadsDescriptorRegisteredIteration:
    """Run the selected order: ``registerShapes`` then control and branching.

    ``atlas_shape_priors`` remain in atlas/template space for registration;
    the transition's ``thinning`` already consumed the registered image-space
    shape values. The same 13-class membership result flows into registration
    and into the control branch. This function composes one iteration only and
    does not update atlas, distance, or correction-field inputs for another.
    """
    if not isinstance(transition, LesionToadsDescriptorContinuingIteration):
        raise TypeError("transition must be a LesionToadsDescriptorContinuingIteration")
    registration = register_after_next_gm_membership_update(
        transition.membership,
        atlas_shape_priors,
        initial_transform,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        align_shapes=align_shapes,
        init_align_iter=init_align_iter,
    )
    thinning_state = transition.thinning.run
    control = continue_after_next_iteration_registration(
        registration,
        segmentation=thinning_state.segmentation,
        skeleton=thinning_state.skeleton,
        ordering=thinning_state.ordering,
        class_volume=thinning_state.class_volume,
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
        retain_growth_event_states=retain_growth_event_states,
    )
    return LesionToadsDescriptorRegisteredIteration(
        transition,
        registration,
        control,
        control.stage_order,
    )
