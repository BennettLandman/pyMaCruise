"""Continue the selected registered TOADS pass through its outer-loop branch.

This composes the P18cyp rigid registration/control handoff with the selected
post-control stop or continue branch. The caller supplies the working-grid
state used by ``AlgorithmMultipleToads``; field refitting and full workflow
assembly remain separate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

from numpy.typing import ArrayLike

from .toads_initial_registration_control import ToadsInitialRegistrationControl
from .toads_iteration import (
    ToadsContinuingTransition,
    ToadsStoppingTransition,
    toads_continuing_iteration,
    toads_stopping_iteration,
)
from .toads_control import toads_iteration_control


ToadsRegisteredIterationBranch: TypeAlias = ToadsStoppingTransition | ToadsContinuingTransition


@dataclass(frozen=True)
class ToadsDefaultRegisteredIterationTransition:
    """Registered membership plus its source-selected outer-loop result."""

    registration_control: ToadsInitialRegistrationControl
    branch: ToadsRegisteredIterationBranch
    stage_order: tuple[str, ...]


def continue_registered_toads_iteration(
    registration_control: ToadsInitialRegistrationControl,
    previous_energy: float,
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    ordering: ArrayLike,
    class_volumes: ArrayLike,
    *,
    do_growing: bool,
    template_labels: ArrayLike | None = None,
    previous_classification: ArrayLike | None = None,
    boundary_distances: ArrayLike | None = None,
    relation_limit: int | None = None,
    available: ArrayLike | None = None,
) -> ToadsDefaultRegisteredIterationTransition:
    """Apply growth/classification or growth/relation refresh after registration.

    The decision is recomputed from the supplied energy and wrapper settings,
    then required to match P18cyp's stored decision. This guards against
    accidentally composing a registration handoff with control parameters
    from a different iteration. The membership map is reused by reference;
    this function neither recomputes it nor transforms its grid.
    """
    if not isinstance(registration_control, ToadsInitialRegistrationControl):
        raise TypeError("registration_control must be a selected P18cyp result")
    membership_pass = registration_control.handoff.membership_pass
    selected_control = toads_iteration_control(
        previous_energy,
        membership_pass.objective,
        iteration,
        maximum_iterations,
        convergence_threshold,
        do_growing=do_growing,
    )
    if selected_control != registration_control.control:
        raise ValueError("outer-loop settings do not reproduce the P18cyp control decision")

    if selected_control.stop:
        if template_labels is None or previous_classification is None:
            raise ValueError("stopping branch requires template_labels and previous_classification")
        branch = toads_stopping_iteration(
            previous_energy,
            membership_pass.objective,
            iteration,
            maximum_iterations,
            convergence_threshold,
            segmentation,
            skeleton,
            membership_pass.memberships,
            ordering,
            class_volumes,
            template_labels,
            previous_classification,
            do_growing=do_growing,
            available=available,
        )
        branch_stages = (
            *(("propagateGrowing",) if selected_control.run_growing else ()),
            "updateClassificationFromLabels",
        )
    else:
        if boundary_distances is None or relation_limit is None:
            raise ValueError("continuing branch requires boundary_distances and relation_limit")
        branch = toads_continuing_iteration(
            previous_energy,
            membership_pass.objective,
            iteration,
            maximum_iterations,
            convergence_threshold,
            segmentation,
            skeleton,
            membership_pass.memberships,
            ordering,
            class_volumes,
            boundary_distances,
            relation_limit=relation_limit,
            do_growing=do_growing,
            available=available,
        )
        branch_stages = (
            "propagateGrowing",
            "computeRelations",
            "computeBoundaryDistances",
            "computeClosestLabels",
        )

    return ToadsDefaultRegisteredIterationTransition(
        registration_control,
        branch,
        (*registration_control.stage_order, *branch_stages),
    )
