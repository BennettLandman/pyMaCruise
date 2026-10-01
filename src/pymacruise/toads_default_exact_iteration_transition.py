"""Connect an exact TOADS membership result to the selected outer-loop branch.

This seam starts after ``computeExactMemberships()`` has produced memberships
and energy. The caller supplies fixed post-registration working-grid state:
the atlas transform, growing ordering, and initial boundary distances. It then
follows ``AlgorithmMultipleToads`` source order through stop/grow/classify or
continue/grow/relation-refresh. Atlas registration itself is outside this
composition and is not validated here.
"""

from __future__ import annotations

from typing import TypeAlias

from numpy.typing import ArrayLike

from .toads_default_atlas_membership import ToadsDefaultAtlasMembershipPass
from .toads_iteration import (
    ToadsContinuingTransition,
    ToadsStoppingTransition,
    toads_continuing_iteration,
    toads_stopping_iteration,
)


ToadsDefaultExactIterationTransition: TypeAlias = (
    ToadsStoppingTransition | ToadsContinuingTransition
)


def toads_default_exact_iteration_transition(
    membership_pass: ToadsDefaultAtlasMembershipPass,
    previous_energy: float,
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    ordering: ArrayLike,
    class_volumes: ArrayLike,
    *,
    template_labels: ArrayLike | None = None,
    previous_classification: ArrayLike | None = None,
    boundary_distances: ArrayLike | None = None,
    relation_limit: int | None = None,
    do_growing: bool,
    available: ArrayLike | None = None,
) -> ToadsDefaultExactIterationTransition:
    """Run the selected stop or continue branch from a fixed exact pass.

    ``membership_pass.objective`` is the current energy; its updated
    memberships seed ``propagateGrowing``. The stop branch transfers labels
    only after its optional final growth. The continuing branch always grows,
    then rebuilds relations, boundary distances, and closest labels. For a
    continuing branch, ``boundary_distances`` represents the fixed working
    state before this iteration's relation refresh.

    Registration is deliberately not called: segmentation, skeleton,
    ordering, and boundary-distance arrays must already share the grid and
    atlas-transform state used by the membership pass. This validates a
    bounded fixed-registration transition, not registration or a complete
    TOADS iteration.
    """
    if not isinstance(membership_pass, ToadsDefaultAtlasMembershipPass):
        raise TypeError("membership_pass must be a selected exact atlas membership result")

    # Match AlgorithmMultipleToads: compute energy ratio and stop decision
    # before deciding which branch runs. Existing branch adapters repeat this
    # source-faithful decision while applying the branch-specific methods.
    from .toads_control import toads_iteration_control

    control = toads_iteration_control(
        previous_energy,
        membership_pass.objective,
        iteration,
        maximum_iterations,
        convergence_threshold,
        do_growing=do_growing,
    )
    if control.stop:
        if template_labels is None or previous_classification is None:
            raise ValueError("stopping branch requires template_labels and previous_classification")
        return toads_stopping_iteration(
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

    if boundary_distances is None:
        raise ValueError("continuing branch requires initial boundary_distances")
    if relation_limit is None:
        raise ValueError("continuing branch requires relation_limit")
    return toads_continuing_iteration(
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
