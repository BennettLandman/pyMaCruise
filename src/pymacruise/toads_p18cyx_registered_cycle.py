"""Continue a P18cyw exact membership update through registration and branch.

This bounded second-loop adapter reuses the selected registration/control and
post-registration branch adapters. Atlas shapes and hard-label state are
explicit fixed synthetic inputs; no new registration, growth, or relation
kernel is implemented here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from .toads_default_membership_registration_control import (
    ToadsDefaultMembershipRegistrationControl,
    continue_default_membership_through_registration_control,
)
from .toads_default_membership_registered_iteration import (
    ToadsDefaultMembershipRegisteredIteration,
    continue_default_membership_registered_iteration,
)
from .toads_p18cyw_exact_iteration import ToadsP18cywExactIteration


@dataclass(frozen=True)
class ToadsP18cyxRegisteredCycle:
    """P18cyw memberships with their next registration and branch state."""

    exact_iteration: ToadsP18cywExactIteration
    registration_control: ToadsDefaultMembershipRegistrationControl
    registered_iteration: ToadsDefaultMembershipRegisteredIteration
    stage_order: tuple[str, ...]


def continue_p18cyw_membership_through_registration_branch(
    exact_iteration: ToadsP18cywExactIteration,
    atlas_shapes: Any,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    do_growing: bool,
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    ordering: ArrayLike,
    class_volumes: ArrayLike,
    boundary_distances: ArrayLike,
    template_labels: ArrayLike,
    previous_classification: ArrayLike,
    relation_limit: int,
    initial_available: ArrayLike | None = None,
    initial_transform: tuple[float, float, float, float, float, float] | None = None,
    align_shapes: bool = True,
    offset: int = 0,
    subsample: int = 2,
    max_registration_iterations: int = 5,
    registration_energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
) -> ToadsP18cyxRegisteredCycle:
    """Run P18cyw energy/registration/stop/grow-or-refresh using existing code.

    The previous objective and rigid transform are taken from the immediately
    preceding P18cyv registered branch, matching the wrapper's carried
    iteration state. The caller supplies the current template atlas and
    post-thinning hard-label/relation state explicitly.
    """
    if not isinstance(exact_iteration, ToadsP18cywExactIteration):
        raise TypeError("exact_iteration must be a selected P18cyw result")
    previous_registration = exact_iteration.p18cyv.registered_iteration.registration_control
    previous_objective = np.float32(
        exact_iteration.p18cyv.registered_iteration.branch.handoff.membership_pass.objective
    )
    transform = (
        tuple(float(value) for value in previous_registration.handoff.transform)
        if initial_transform is None
        else initial_transform
    )
    current_membership = exact_iteration.exact_iteration.membership_pass

    registration_control = continue_default_membership_through_registration_control(
        current_membership,
        atlas_shapes,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        initial_transform=transform,
        previous_energy=previous_objective,
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
        previous_energy=previous_objective,
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
    stages = (
        *exact_iteration.stage_order,
        *registration_control.stage_order[1:],
        *registered_iteration.stage_order[len(registration_control.stage_order):],
    )
    return ToadsP18cyxRegisteredCycle(
        exact_iteration=exact_iteration,
        registration_control=registration_control,
        registered_iteration=registered_iteration,
        stage_order=tuple(stages),
    )
