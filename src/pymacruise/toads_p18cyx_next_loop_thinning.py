"""Continue a P18cyx continuing branch into the next thinning pass.

This reuses the selected P18cyv short-constructor thinning adapter after the
second-loop registration and growth/relation refresh. It introduces no
thinning or topology kernel.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from .toads_post_registration_next_thinning import (
    ToadsPostRegistrationNextThinning,
    toads_post_registration_next_thinning,
)
from .toads_p18cyx_registered_cycle import ToadsP18cyxRegisteredCycle


@dataclass(frozen=True)
class ToadsP18cyxNextLoopThinning:
    """P18cyx branch state advanced through the next short-constructor pass."""

    registered_cycle: ToadsP18cyxRegisteredCycle
    next_thinning: ToadsPostRegistrationNextThinning
    stage_order: tuple[str, ...]


def continue_p18cyx_branch_through_next_thinning(
    registered_cycle: ToadsP18cyxRegisteredCycle,
    transformed_shape_values: ArrayLike,
) -> ToadsP18cyxNextLoopThinning:
    """Run the selected next-loop thinning kernel on the current branch output."""
    if not isinstance(registered_cycle, ToadsP18cyxRegisteredCycle):
        raise TypeError("registered_cycle must be a selected P18cyx result")
    branch = registered_cycle.registered_iteration.branch
    shapes = np.asarray(transformed_shape_values)
    memberships = np.asarray(branch.handoff.membership_pass.memberships)
    if shapes.shape != memberships.shape:
        raise ValueError("transformed_shape_values must match current membership maps")
    result = toads_post_registration_next_thinning(branch, shapes)
    return ToadsP18cyxNextLoopThinning(
        registered_cycle,
        result,
        (*registered_cycle.stage_order, "next_iteration.propagateThinning"),
    )
