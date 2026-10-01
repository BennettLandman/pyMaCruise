"""Continue a bounded P18cyt registered iteration through next-loop thinning.

This composes the existing P18cyu branch with the selected short-constructor
thinning adapter. Membership correction, transformed shapes, and registration
remain fixed inputs at this seam; it is not a complete TOADS loop.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from .toads_default_membership_registered_iteration import (
    ToadsDefaultMembershipRegisteredIteration,
)
from .toads_post_registration_next_thinning import (
    ToadsPostRegistrationNextThinning,
    toads_post_registration_next_thinning,
)


@dataclass(frozen=True)
class ToadsP18cyvNextIteration:
    """Registered iteration and its selected next thinning pass."""

    registered_iteration: ToadsDefaultMembershipRegisteredIteration
    next_thinning: ToadsPostRegistrationNextThinning
    stage_order: tuple[str, ...]


def continue_p18cyv_iteration_through_next_thinning(
    registered_iteration: ToadsDefaultMembershipRegisteredIteration,
    transformed_shape_values: ArrayLike,
) -> ToadsP18cyvNextIteration:
    """Run only the selected next-loop exact-membership thinning pass."""
    if not isinstance(registered_iteration, ToadsDefaultMembershipRegisteredIteration):
        raise TypeError("registered_iteration must be a selected P18cyu continuation")
    shapes = np.asarray(transformed_shape_values)
    memberships = registered_iteration.branch.handoff.membership_pass.memberships
    if shapes.shape != memberships.shape:
        raise ValueError("transformed_shape_values must match the current membership maps")

    result = toads_post_registration_next_thinning(
        registered_iteration.branch,
        shapes,
    )
    if result.branch is not registered_iteration.branch:
        raise RuntimeError("next-loop thinning changed the registered branch handoff")
    stages = (*registered_iteration.stage_order, "next_iteration.propagateThinning")
    return ToadsP18cyvNextIteration(
        registered_iteration,
        result,
        stages,
    )
