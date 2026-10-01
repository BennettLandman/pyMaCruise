"""Continue P18czf's non-stopping branch through the next thinning pass.

This named P18 continuation delegates to the existing selected
post-registration thinning kernel with current memberships and refreshed
relation adjacency. It introduces no new thinning or topology logic.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from .toads_p18czf_registered_branch import ToadsP18czfRegisteredBranch
from .toads_post_registration_next_thinning import (
    ToadsPostRegistrationNextThinning,
    toads_post_registration_next_thinning,
)


@dataclass(frozen=True)
class ToadsP18czgNextLoopThinning:
    """P18czf branch state advanced through the next short-constructor pass."""

    registered_branch: ToadsP18czfRegisteredBranch
    next_thinning: ToadsPostRegistrationNextThinning
    stage_order: tuple[str, ...]


def continue_p18czf_branch_through_next_loop_thinning(
    registered_branch: ToadsP18czfRegisteredBranch,
    transformed_shape_values: ArrayLike,
) -> ToadsP18czgNextLoopThinning:
    """Run one exact-membership ``propagateThinning`` after P18czf."""
    if not isinstance(registered_branch, ToadsP18czfRegisteredBranch):
        raise TypeError("registered_branch must be a selected P18czf result")
    branch = registered_branch.registered_iteration.branch
    shapes = np.asarray(transformed_shape_values)
    memberships = np.asarray(branch.handoff.membership_pass.memberships)
    if shapes.shape != memberships.shape:
        raise ValueError("transformed_shape_values must match the registered memberships")
    thinning = toads_post_registration_next_thinning(branch, shapes)
    return ToadsP18czgNextLoopThinning(
        registered_branch,
        thinning,
        (*registered_branch.stage_order, "next_iteration.propagateThinning"),
    )
