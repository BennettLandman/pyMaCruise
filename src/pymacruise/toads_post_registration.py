"""Selected outer-loop decision and state handoff after TOADS registration.

This composes already validated fixed-input kernels. Membership computation
and atlas registration are supplied as a ``ToadsDefaultRegistrationHandoff``;
this function does not repeat either stage or introduce another iteration.
Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``. This
captures only the selected wrapper's post-registration control-flow branch,
not the complete topology-preserving evolution described in the papers.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_control import (
    ToadsIterationControl,
    toads_classification_from_labels,
    toads_iteration_control,
)
from .toads_default_registration_handoff import ToadsDefaultRegistrationHandoff
from .toads_relations import ToadsRelationState, toads_refresh_relation_state
from .toads_thinning import ToadsGrowingPass, toads_propagate_growing


@dataclass(frozen=True)
class ToadsPostRegistrationBranch:
    """Result of one source-ordered post-registration wrapper branch."""

    handoff: ToadsDefaultRegistrationHandoff
    control: ToadsIterationControl
    growing: ToadsGrowingPass | None
    classification: NDArray[np.int8] | None
    changed_classification_voxels: int | None
    relations: ToadsRelationState | None
    stage_order: tuple[str, ...]


def toads_post_registration_branch(
    handoff: ToadsDefaultRegistrationHandoff,
    *,
    previous_energy: float,
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    ordering: ArrayLike,
    class_volumes: ArrayLike,
    boundary_distances: ArrayLike,
    template_labels: ArrayLike,
    previous_classification: ArrayLike,
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    relation_limit: int,
    do_growing: bool,
    initial_available: ArrayLike | None = None,
    record_growing_events: bool = True,
) -> ToadsPostRegistrationBranch:
    """Run only the selected wrapper branch immediately after registration.

    The current energy and 11-class memberships come from ``handoff``. The
    source compares energy, checks convergence and the iteration ceiling,
    grows when ``(!stop) || doGrowing``, then either transfers the final hard
    labels or refreshes relations, boundary distances, and closest labels.
    The selected Java wrapper performs registration before these decisions.
    """
    if not isinstance(handoff, ToadsDefaultRegistrationHandoff):
        raise TypeError("handoff must be a ToadsDefaultRegistrationHandoff")

    segment = np.asarray(segmentation)
    mem = np.asarray(handoff.membership_pass.memberships)
    if segment.ndim != 3 or mem.ndim != 4 or mem.shape[:3] != segment.shape:
        raise ValueError("segmentation and registered membership maps must share a 3-D grid")
    if mem.shape[-1] != 11:
        raise ValueError("the post-registration default-atlas branch requires 11 membership classes")

    stages = ["membership_energy", "registration_complete", "energy_ratio", "stop_checks"]
    control = toads_iteration_control(
        previous_energy,
        handoff.membership_pass.objective,
        iteration,
        maximum_iterations,
        convergence_threshold,
        do_growing=do_growing,
    )

    growing = None
    if control.run_growing:
        stages.append("propagateGrowing")
        growing = toads_propagate_growing(
            segment,
            skeleton,
            mem,
            ordering,
            class_volumes,
            initial_available=initial_available,
            record_events=record_growing_events,
        )
        segment_after_growth = growing.segmentation
    else:
        segment_after_growth = segment

    classification = None
    changed = None
    relations = None
    if control.update_classification:
        stages.append("updateClassificationFromLabels")
        classification, changed = toads_classification_from_labels(
            segment_after_growth,
            template_labels,
            previous_classification,
        )
    elif control.refresh_relations:
        stages.extend(("computeRelations", "computeBoundaryDistances", "computeClosestLabels"))
        relations = toads_refresh_relation_state(
            segment_after_growth,
            mem,
            boundary_distances,
            relation_limit,
        )

    return ToadsPostRegistrationBranch(
        handoff=handoff,
        control=control,
        growing=growing,
        classification=classification,
        changed_classification_voxels=changed,
        relations=relations,
        stage_order=tuple(stages),
    )
