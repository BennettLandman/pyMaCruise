"""Continue a non-stopping registered TOADS branch into the next loop pass.

The selected ``AlgorithmMultipleToads`` wrapper refreshes relation state on a
continuing iteration, returns to the loop, and starts the next iteration with
``propagateThinning()``. This adapter passes the updated labels, current exact
memberships, and refreshed class adjacency to the selected short-constructor
thinning defaults. Registration, transformed shape generation, and any
inhomogeneity-field refit remain fixed inputs outside this slice.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from .toads_post_registration import ToadsPostRegistrationBranch
from .toads_thinning import ToadsThinningPass, toads_propagate_thinning


@dataclass(frozen=True)
class ToadsPostRegistrationNextThinning:
    """The source branch and its next-loop thinning result."""

    branch: ToadsPostRegistrationBranch
    thinning: ToadsThinningPass
    stage_order: tuple[str, ...]


def toads_post_registration_next_thinning(
    branch: ToadsPostRegistrationBranch,
    transformed_shape_values: ArrayLike,
    *,
    record_thinning_items: bool = True,
) -> ToadsPostRegistrationNextThinning:
    """Run one next-loop ``propagateThinning`` after a continuing branch.

    The selected ``MedicAlgorithmToads08`` short constructor sets
    ``useApprox=false``; the delegated thinning kernel therefore uses exact
    class memberships. It also receives the selected defaults
    ``volumeRatio=0.5`` and ``limit=3``. Correction fields and transformed
    shapes are held fixed at this handoff.
    """
    if not isinstance(branch, ToadsPostRegistrationBranch):
        raise TypeError("branch must be a selected post-registration result")
    if branch.control.stop or branch.growing is None or branch.relations is None:
        raise ValueError("next-loop thinning requires a non-stopping branch with growth and refreshed relations")

    memberships = np.asarray(branch.handoff.membership_pass.memberships)
    labels = np.asarray(branch.growing.segmentation)
    shapes = np.asarray(transformed_shape_values)
    adjacency = branch.relations.relation.class_adjacency
    if memberships.ndim != 4 or memberships.shape[:3] != labels.shape:
        raise ValueError("registered memberships and grown labels must share a 3-D grid")
    if shapes.shape != memberships.shape:
        raise ValueError("transformed_shape_values must match the registered membership maps")
    if adjacency.shape != (memberships.shape[-1], memberships.shape[-1]):
        raise ValueError("refreshed class adjacency must match the membership class count")

    thinning = toads_propagate_thinning(
        labels,
        memberships,
        adjacency,
        shapes,
        volume_ratio=0.5,
        limit=3.0,
        record_processed_items=record_thinning_items,
    )
    return ToadsPostRegistrationNextThinning(
        branch=branch,
        thinning=thinning,
        stage_order=(
            *branch.stage_order,
            "next_iteration.propagateThinning",
        ),
    )
