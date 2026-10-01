"""Continue P18czc's non-stopping branch through the next thinning pass.

This labeled adapter delegates the source-ordered pass to the existing
post-registration thinning kernel. It carries current exact memberships,
grown labels, and refreshed class adjacency together; no new thinning or
topology algorithm is introduced here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from .toads_p18czc_registered_branch import ToadsP18czcRegisteredBranch
from .toads_post_registration_next_thinning import (
    ToadsPostRegistrationNextThinning,
    toads_post_registration_next_thinning,
)


@dataclass(frozen=True)
class ToadsP18czdNextLoopThinning:
    """P18czc registered branch and its next short-constructor thinning pass."""

    registered_branch: ToadsP18czcRegisteredBranch
    next_thinning: ToadsPostRegistrationNextThinning
    stage_order: tuple[str, ...]


def continue_p18czc_branch_through_next_loop_thinning(
    registered_branch: ToadsP18czcRegisteredBranch,
    transformed_shape_values: ArrayLike,
) -> ToadsP18czdNextLoopThinning:
    """Run one exact-membership thinning pass after P18czc relation refresh."""
    if not isinstance(registered_branch, ToadsP18czcRegisteredBranch):
        raise TypeError("registered_branch must be a selected P18czc registered branch")
    branch = registered_branch.registered_iteration.branch
    shapes = np.asarray(transformed_shape_values)
    memberships = np.asarray(branch.handoff.membership_pass.memberships)
    if shapes.shape != memberships.shape:
        raise ValueError("transformed_shape_values must match the registered memberships")
    result = toads_post_registration_next_thinning(branch, shapes)
    return ToadsP18czdNextLoopThinning(
        registered_branch,
        result,
        (*registered_branch.stage_order, "next_iteration.propagateThinning"),
    )
