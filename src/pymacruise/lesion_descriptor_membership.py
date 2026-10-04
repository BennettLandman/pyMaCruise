"""Connect descriptor topology labels to the selected membership updater.

In the active VandyCRUISE membership loop, OBJECT, MASK, and active
OPTIMIZED_CSF (whose data term falls through to ``outlier``) use the ordinary
class smoothing, prior, inversion, and normalization update. Their distinct
raw data formulas are calculated upstream. LESION additionally updates the
shared lesion-membership map. This module preserves that split between class
topology and already-computed data terms.
"""

from __future__ import annotations

from typing import Sequence

from .lesion_objective_membership_pass import (
    LesionObjectiveMembershipPass,
    run_lesion_toads_objective_membership_pass,
)


_MEMBERSHIP_TYPE = {
    "obj": "OBJECT",
    "object": "OBJECT",
    "mask": "MASK",
    "optimized_csf": "MASK",
    "lesion": "LESION",
}


def map_lesion_descriptor_membership_types(topologies: Sequence[str]) -> tuple[str, ...]:
    """Map the bundled descriptor families to active update categories.

    ``optimized_csf`` maps to MASK only for membership-update dispatch: the
    selected source does not run the commented optimized-CSF data formula,
    and both categories take the ordinary class membership branch. Its raw
    energy must still be produced from the source fallback before calling the
    membership pass. Unrecognized families are rejected rather than folded
    into a generic approximation.
    """
    if not topologies or any(not isinstance(value, str) for value in topologies):
        raise ValueError("topologies must be a nonempty sequence of strings")
    mapped: list[str] = []
    for value in topologies:
        token = value.lower()
        if token not in _MEMBERSHIP_TYPE:
            raise ValueError(f"unsupported LesionTOADS descriptor membership topology: {value}")
        mapped.append(_MEMBERSHIP_TYPE[token])
    if mapped.count("LESION") != 1:
        raise ValueError("descriptor membership requires exactly one LESION class")
    return tuple(mapped)


def run_lesion_toads_descriptor_membership_pass(
    *args: object,
    descriptor_types: Sequence[str],
    **kwargs: object,
) -> LesionObjectiveMembershipPass:
    """Run the in-place membership kernel with descriptor-derived categories.

    All array/term arguments are the fixed inputs accepted by
    :func:`run_lesion_toads_objective_membership_pass`; ``descriptor_types``
    replaces its lower-level ``class_types`` argument.
    """
    if "class_types" in kwargs:
        raise TypeError("pass descriptor_types instead of class_types")
    mapped_types = map_lesion_descriptor_membership_types(descriptor_types)
    return run_lesion_toads_objective_membership_pass(
        *args,
        class_types=mapped_types,
        **kwargs,  # type: ignore[arg-type]
    )
