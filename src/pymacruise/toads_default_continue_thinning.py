"""Compose the selected continuing TOADS branch into the next thinning pass.

This transition starts with P18cyj's fixed-registration continuation state.
It takes the updated labels and exact membership map, uses the newly refreshed
class adjacency, and invokes the selected wrapper's next-iteration thinning
defaults. Registration remains fixed and outside this composition.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from .toads_iteration import ToadsContinuingTransition
from .toads_thinning import ToadsThinningPass, toads_propagate_thinning


def toads_default_continue_thinning(
    transition: ToadsContinuingTransition,
    memberships: ArrayLike,
    transformed_shape_values: ArrayLike,
) -> ToadsThinningPass:
    """Run ``initThinningLabels`` and ``propagateThinning`` after continuation.

    ``AlgorithmMultipleToads`` returns to the top of its loop after a
    non-stopping branch has grown labels and refreshed relation state; the
    first next-iteration operation is ``propagateThinning``. The selected
    short ``MedicAlgorithmToads08`` constructor uses exact memberships,
    ``volumeRatio=0.5``, and ``limit=3``. Inputs representing registration and
    transformed shape values are fixed at this seam; no registration call is
    made here.
    """
    if not isinstance(transition, ToadsContinuingTransition):
        raise TypeError("transition must be a selected continuing-iteration result")
    labels = transition.growing.segmentation
    relation = transition.relation_state.relation
    mems = np.asarray(memberships)
    shapes = np.asarray(transformed_shape_values)
    if mems.ndim != 4 or tuple(mems.shape[:3]) != tuple(labels.shape):
        raise ValueError("memberships must match the continued segmentation grid")
    if mems.shape[-1] != relation.class_count:
        raise ValueError("memberships and refreshed relation state must have the same class count")
    if shapes.ndim != 4 or tuple(shapes.shape) != tuple(mems.shape):
        raise ValueError("transformed_shape_values must match memberships")
    adjacency = relation.class_adjacency
    if adjacency.shape != (relation.class_count, relation.class_count):
        raise ValueError("refreshed relation state does not provide a class adjacency matrix")
    return toads_propagate_thinning(
        labels,
        mems,
        adjacency,
        shapes,
        volume_ratio=0.5,
        limit=3.0,
    )
