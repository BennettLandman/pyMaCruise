"""Continue P18czi's registered branch into the next thinning pass.

The adapter reuses the exact-membership short-constructor thinning kernel with
the current registered memberships, newly grown labels, and refreshed
adjacency. Shape values remain explicit synthetic inputs; no thinning kernel
is added here.
"""

from __future__ import annotations

from dataclasses import dataclass

from numpy.typing import ArrayLike

from .toads_p18czi_registered_branch import ToadsP18cziRegisteredBranch
from .toads_post_registration_next_thinning import (
    ToadsPostRegistrationNextThinning,
    toads_post_registration_next_thinning,
)


@dataclass(frozen=True)
class ToadsP18czjNextLoopThinning:
    """P18czi registered branch plus its following selected thinning result."""

    registered_branch: ToadsP18cziRegisteredBranch
    next_thinning: ToadsPostRegistrationNextThinning
    stage_order: tuple[str, ...]


def continue_p18czi_branch_through_next_loop_thinning(
    registered_branch: ToadsP18cziRegisteredBranch,
    transformed_shape_values: ArrayLike,
) -> ToadsP18czjNextLoopThinning:
    """Run one selected short-constructor thinning pass after P18czi."""
    if not isinstance(registered_branch, ToadsP18cziRegisteredBranch):
        raise TypeError("registered_branch must be a selected P18czi result")
    branch = registered_branch.registered_iteration.branch
    thinning = toads_post_registration_next_thinning(branch, transformed_shape_values)
    return ToadsP18czjNextLoopThinning(
        registered_branch,
        thinning,
        (*registered_branch.stage_order, "next_iteration.propagateThinning"),
    )
