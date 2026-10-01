"""Carry P18czk memberships through the configured final registration/control.

The selected rigid registration and wrapper stop/classification decision are
reused unchanged. The returned branch is intentionally terminal for this
synthetic loop because the selected iteration count reaches nIterMax.
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
from .toads_p18czk_next_exact_iteration import ToadsP18czkNextExactIteration


@dataclass(frozen=True)
class ToadsP18czlTerminalBranch:
    """P18czk memberships after final selected registration/control handling."""

    exact_iteration: ToadsP18czkNextExactIteration
    previous_energy: np.float32
    initial_transform: np.ndarray
    registration_control: ToadsDefaultMembershipRegistrationControl
    registered_iteration: ToadsDefaultMembershipRegisteredIteration
    stage_order: tuple[str, ...]


def continue_p18czk_membership_through_terminal_branch(
    exact_iteration: ToadsP18czkNextExactIteration,
    atlas_shapes: Any,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    template_labels: ArrayLike,
    previous_classification: ArrayLike,
    relation_limit: int,
    align_shapes: bool = True,
    offset: int = 0,
    subsample: int = 2,
    max_registration_iterations: int = 5,
    registration_energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
) -> ToadsP18czlTerminalBranch:
    """Run final registration/control with growth disabled by the iteration cap."""
    if not isinstance(exact_iteration, ToadsP18czkNextExactIteration):
        raise TypeError("exact_iteration must be a selected P18czk result")
    previous_loop = exact_iteration.previous_loop
    prior_registered = previous_loop.registered_branch
    prior_branch = previous_loop.next_thinning.branch
    previous_energy = np.float32(prior_branch.handoff.membership_pass.objective)
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
        do_growing=False,
        align_shapes=align_shapes,
        offset=offset,
        subsample=subsample,
        max_registration_iterations=max_registration_iterations,
        registration_energy_tolerance=registration_energy_tolerance,
        minimum_lambda=minimum_lambda,
    )

    thinning = previous_loop.next_thinning.thinning
    branch = previous_loop.next_thinning.branch
    if branch.growing is None or branch.relations is None:
        raise ValueError("terminal branch requires carried growth and refreshed relation state")
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
        do_growing=False,
        initial_available=thinning.available,
    )
    if not registered_iteration.branch.control.stop:
        raise RuntimeError("configured final iteration did not produce a terminal wrapper decision")
    return ToadsP18czlTerminalBranch(
        exact_iteration,
        previous_energy,
        initial_transform,
        registration_control,
        registered_iteration,
        registered_iteration.stage_order,
    )
