"""State handoff for the selected LesionTOADS transformed topology.

This mirrors ``LesionToadPreprocess.updateTransformedTemplate``: replace the
preprocessing topology reference with the just-resampled atlas and set its
``transformed`` flag. The selected wrapper invokes it immediately after
``LesionToadAtlas.computeTransformedTemplate`` on the optional
``alignShapes`` path. The input reference is retained without copying, as in
the Java assignment.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray


@dataclass
class LesionToadsTemplateState:
    """Only the two mutable fields affected by the Java method."""

    topology: NDArray[Any] | None = None
    transformed: bool = False


def update_lesion_transformed_template(
    state: LesionToadsTemplateState,
    transformed_template: NDArray[Any],
) -> None:
    """Assign the transformed template reference and mark it active.

    This intentionally performs no copy, resampling, label validation, or
    geometry change; the source method consists only of the two assignments
    ``topology = tpl`` and ``transformed = true``.
    """
    state.topology = transformed_template
    state.transformed = True
