"""Continue P18cze memberships through the next registration/control branch.

The previous energy and rigid transform are carried from the prior P18czb
registration. The current memberships/objective come from P18cze; thinning and
relation inputs are carried from P18czd/P18czc. This is an explicit synthetic
single-channel continuation over existing adapters, not a complete workflow.
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
from .toads_p18cze_next_exact_iteration import ToadsP18czeNextExactIteration


@dataclass(frozen=True)
class ToadsP18czfRegisteredBranch:
    """P18cze exact memberships after carried registration and loop control."""

    exact_iteration: ToadsP18czeNextExactIteration
    previous_energy: np.float32
    initial_transform: np.ndarray
    registration_control: ToadsDefaultMembershipRegistrationControl
    registered_iteration: ToadsDefaultMembershipRegisteredIteration
    stage_order: tuple[str, ...]


def continue_p18cze_membership_through_registered_branch(
    exact_iteration: ToadsP18czeNextExactIteration,
    atlas_shapes: Any,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    do_growing: bool,
    template_labels: ArrayLike,
    previous_classification: ArrayLike,
    relation_limit: int,
    align_shapes: bool = True,
    offset: int = 0,
    subsample: int = 2,
    max_registration_iterations: int = 5,
    registration_energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
) -> ToadsP18czfRegisteredBranch:
    """Carry the prior objective/transform and current next-loop state forward.

    The P18cza exact objective is the previous energy for P18cze. The transform
    is the registered transform from the preceding P18czb transition. P18czd
    supplies the just-thinned labels, skeleton, availability, and volumes;
    P18czc supplies the carried ordering and boundary distances.
    """
    if not isinstance(exact_iteration, ToadsP18czeNextExactIteration):
        raise TypeError("exact_iteration must be a selected P18cze result")
    previous_thinning = exact_iteration.previous_thinning
    prior_registered = previous_thinning.registered_branch
    prior_exact = prior_registered.exact_iteration.exact_iteration
    previous_energy = np.float32(prior_exact.membership_pass.objective)
    initial_transform = np.asarray(
        prior_registered.registration_control.handoff.transform, dtype=np.float32
    ).copy()
    current = exact_iteration.exact_iteration.exact_iteration.membership_pass

    registration_control = continue_default_membership_through_registration_control(
        current,
        atlas_shapes,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        initial_transform=tuple(float(v) for v in initial_transform),
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
    thinning = previous_thinning.next_thinning.thinning
    branch = previous_thinning.next_thinning.branch
    if branch.growing is None or branch.relations is None:
        raise ValueError("P18czf requires P18czd continuing growth and refreshed relation state")
    registered_iteration = continue_default_membership_registered_iteration(
        registration_control,
        previous_energy=previous_energy,
        iteration=iteration,
        maximum_iterations=maximum_iterations,
        convergence_threshold=convergence_threshold,
        segmentation=thinning.segmentation,
        skeleton=thinning.skeleton,
        ordering=branch.growing.ordering,
        class_volumes=thinning.class_volumes,
        boundary_distances=branch.relations.boundary_distances,
        template_labels=template_labels,
        previous_classification=previous_classification,
        relation_limit=relation_limit,
        do_growing=do_growing,
        initial_available=thinning.available,
    )
    stages = (
        *exact_iteration.stage_order,
        *registration_control.stage_order[1:],
        *registered_iteration.stage_order[len(registration_control.stage_order):],
    )
    return ToadsP18czfRegisteredBranch(
        exact_iteration=exact_iteration,
        previous_energy=previous_energy,
        initial_transform=initial_transform,
        registration_control=registration_control,
        registered_iteration=registered_iteration,
        stage_order=tuple(stages),
    )
