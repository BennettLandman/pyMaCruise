"""Initialize relation, boundary-distance, and closest-label TOADS state.

The selected ``AlgorithmMultipleToads`` startup path calls the three
``MultipleToadSegmentation`` preparation methods in order. The segmentation
constructor initializes its class-wise boundary-distance volume to zero; this
helper supplies that initial state and composes the already source-checked
relation refresh operations. It does not compute class variances, transformed
atlas shapes, or the membership update that follows.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from .toads_relations import ToadsRelationState, toads_refresh_relation_state


def toads_initialize_relation_state(
    segmentation: ArrayLike,
    memberships: Any,
    limit: int,
) -> ToadsRelationState:
    """Run startup relation, boundary-distance, and closest-label preparation.

    The selected constructor allocates ``boundaryDist`` as an all-zero
    ``(x,y,z,class)`` signed-byte array before ``AlgorithmMultipleToads`` calls
    ``computeRelations``, ``computeBoundaryDistances``, then
    ``computeClosestLabels``. This wrapper recreates that initial array and
    delegates to the ordered source-matched composition. Input memberships
    and segmentation are the already initialized fixed state; this function
    does not estimate them or implement approximate-membership mode.
    """
    labels = np.asarray(segmentation)
    mems = np.asarray(memberships)
    if labels.ndim != 3 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be a 3-D integer array")
    if mems.ndim != 4 or mems.shape[:3] != labels.shape or mems.shape[-1] < 1:
        raise ValueError("memberships must match the segmentation grid and have a class axis")
    if any(size < 3 for size in labels.shape):
        raise ValueError("segmentation dimensions must be at least three voxels")
    if not np.issubdtype(mems.dtype, np.number) or not np.isfinite(mems).all():
        raise ValueError("memberships must contain finite numeric values")

    initial_boundary_distances = np.zeros(mems.shape, dtype=np.int8)
    return toads_refresh_relation_state(
        labels,
        mems,
        initial_boundary_distances,
        limit,
    )
