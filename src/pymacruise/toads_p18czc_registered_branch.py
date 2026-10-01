"""Continue a bounded P18cza result through registration and loop branching.

This is an orchestration adapter over the existing membership registration,
control, growth, and relation-refresh kernels. All post-thinning hard-label
state is explicit input; it is one synthetic single-channel loop slice, not a
complete TOADS implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from numpy.typing import ArrayLike

from .toads_default_membership_registration_control import (
    ToadsDefaultMembershipRegistrationControl,
    continue_default_membership_through_registration_control,
)
from .toads_default_membership_registered_iteration import (
    ToadsDefaultMembershipRegisteredIteration,
    continue_default_membership_registered_iteration,
)
from .toads_p18cza_carried_degree_iteration import ToadsP18czaCarriedDegreeIteration


@dataclass(frozen=True)
class ToadsP18czcRegisteredBranch:
    """P18cza membership with its registration, wrapper branch, and state order."""

    exact_iteration: ToadsP18czaCarriedDegreeIteration
    registration_control: ToadsDefaultMembershipRegistrationControl
    registered_iteration: ToadsDefaultMembershipRegisteredIteration
    stage_order: tuple[str, ...]


def continue_p18cza_membership_through_registered_branch(
    exact_iteration: ToadsP18czaCarriedDegreeIteration,
    atlas_shapes: Any,
    *,
    previous_energy: float,
    initial_transform: tuple[float, float, float, float, float, float],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    ordering: ArrayLike,
    class_volumes: ArrayLike,
    boundary_distances: ArrayLike,
    template_labels: ArrayLike,
    previous_classification: ArrayLike,
    relation_limit: int,
    do_growing: bool,
    initial_available: ArrayLike | None = None,
    align_shapes: bool = True,
    offset: int = 0,
    subsample: int = 2,
    max_registration_iterations: int = 5,
    registration_energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
) -> ToadsP18czcRegisteredBranch:
    """Run the existing source-ordered registration and post-registration branch.

    The caller supplies the previous loop's membership objective, transform,
    and spatial continuation state. Registration/control retains the exact
    P18cza membership object; the registered iteration verifies the same stop
    and growth decision before invoking the selected growth/relation helpers.
    """
    if not isinstance(exact_iteration, ToadsP18czaCarriedDegreeIteration):
        raise TypeError("exact_iteration must be a bounded P18cza result")
    membership = exact_iteration.exact_iteration.membership_pass
    registration_control = continue_default_membership_through_registration_control(
        membership,
        atlas_shapes,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        initial_transform=initial_transform,
        previous_energy=previous_energy,
        iteration=iteration,
        maximum_iterations=maximum_iterations,
        convergence_threshold=convergence_threshold,
        do_growing=do_growing,
        align_shapes=align_shapes,
        offset=offset,
        subsample=subsample,
        max_registration_iterations=max_registration_iterations,
        registration_energy_tolerance=registration_energy_tolerance,
        minimum_lambda=minimum_lambda,
    )
    registered_iteration = continue_default_membership_registered_iteration(
        registration_control,
        previous_energy=previous_energy,
        iteration=iteration,
        maximum_iterations=maximum_iterations,
        convergence_threshold=convergence_threshold,
        segmentation=segmentation,
        skeleton=skeleton,
        ordering=ordering,
        class_volumes=class_volumes,
        boundary_distances=boundary_distances,
        template_labels=template_labels,
        previous_classification=previous_classification,
        relation_limit=relation_limit,
        do_growing=do_growing,
        initial_available=initial_available,
    )
    return ToadsP18czcRegisteredBranch(
        exact_iteration=exact_iteration,
        registration_control=registration_control,
        registered_iteration=registered_iteration,
        stage_order=(
            *exact_iteration.stage_order,
            *registration_control.stage_order[1:],
            *registered_iteration.stage_order[len(registration_control.stage_order):],
        ),
    )
