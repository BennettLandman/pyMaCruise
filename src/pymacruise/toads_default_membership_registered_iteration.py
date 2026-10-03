"""Continue P18cyt membership/registration control through the selected branch.

This composes the selected wrapper stop/grow/relation decisions with the
existing growth and relation kernels. It adds no image-processing kernel and
does not implement a complete TOADS iteration. Paper context: Bazin and Pham
(2007, 2008), ``docs/papers/toads.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

from numpy.typing import ArrayLike

from .toads_control import toads_iteration_control
from .toads_default_membership_registration_control import (
    ToadsDefaultMembershipRegistrationControl,
)
from .toads_post_registration import (
    ToadsPostRegistrationBranch,
    toads_post_registration_branch,
)


@dataclass(frozen=True)
class ToadsDefaultMembershipRegisteredIteration:
    """P18cyt registered membership state plus its selected loop branch."""

    registration_control: ToadsDefaultMembershipRegistrationControl
    branch: ToadsPostRegistrationBranch
    stage_order: tuple[str, ...]


def continue_default_membership_registered_iteration(
    registration_control: ToadsDefaultMembershipRegistrationControl,
    *,
    previous_energy: float,
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
    record_growing_events: bool = True,
) -> ToadsDefaultMembershipRegisteredIteration:
    """Compose the existing branch adapter after a completed P18cyt slice.

    The existing branch adapter re-evaluates wrapper control from the
    membership objective. Check that decision against the P18cyt registration
    result before accepting its growth/classification/relation state. The
    function holds no extra statistical, topology, growth, or relation logic.
    """
    if not isinstance(registration_control, ToadsDefaultMembershipRegistrationControl):
        raise TypeError("registration_control must be a completed selected P18cyt result")
    handoff = registration_control.handoff
    selected = toads_iteration_control(
        previous_energy,
        handoff.membership_pass.objective,
        iteration,
        maximum_iterations,
        convergence_threshold,
        do_growing=do_growing,
    )
    if selected != registration_control.control:
        raise ValueError("outer-loop settings do not reproduce the P18cyt control decision")

    branch = toads_post_registration_branch(
        handoff,
        previous_energy=previous_energy,
        segmentation=segmentation,
        skeleton=skeleton,
        ordering=ordering,
        class_volumes=class_volumes,
        boundary_distances=boundary_distances,
        template_labels=template_labels,
        previous_classification=previous_classification,
        iteration=iteration,
        maximum_iterations=maximum_iterations,
        convergence_threshold=convergence_threshold,
        relation_limit=relation_limit,
        do_growing=do_growing,
        initial_available=initial_available,
        record_growing_events=record_growing_events,
    )
    if branch.control != registration_control.control:
        raise RuntimeError("post-registration branch changed the selected P18cyt control decision")

    stages = list(registration_control.stage_order)
    if branch.growing is not None:
        stages.append("propagateGrowing")
    if branch.control.update_classification:
        stages.append("updateClassificationFromLabels")
    else:
        stages.extend(("computeRelations", "computeBoundaryDistances", "computeClosestLabels"))
    return ToadsDefaultMembershipRegisteredIteration(
        registration_control,
        branch,
        tuple(stages),
    )
