"""Source-faithful distance support used by LesionTOADS.

This module ports only the small capped distance-to-mask primitive from
``LesionToadSegmentation.fastDistance(boolean[][][], byte, int, int, int)``.
It is not a LesionTOADS implementation or evidence of segmentation
equivalence. The source uses six-connected voxel steps, not physical distances.

Paper context: topology-preserving segmentation is a central part of the
method described by Shiee et al. (2010), but this helper is only one selected
VandyCRUISE distance-support routine. It does not implement the paper's
topological atlas, fast-marching evolution, or lesion segmentation. See
``docs/papers/lesion-toads.md`` for this distinction.

Arrays use ``(x, y, z)`` index order, matching the Java source's
``boolean[x][y][z]`` layout. Callers are responsible for providing arrays in
the source algorithm's RAS-oriented voxel order.
"""

from __future__ import annotations

import operator

import numpy as np


def lesion_toads_distance_to_mask(
    object_mask: np.ndarray, max_distance: int
) -> np.ndarray:
    """Return the selected LesionTOADS six-neighbor capped distance field.

    Object voxels are 0. Non-object voxels have their shortest six-connected
    distance to the object, in voxel steps, through ``max_distance`` inclusive
    for caps of 2 or greater. The source method always executes one propagation
    scan, so a cap of 1 also assigns distance 2 to voxels two steps from the
    object. Voxels beyond the resulting cap and all voxels for an empty object
    remain -1. The result is ``int8``, preserving Java signed-byte output for
    valid source parameter values.

    The selected Java method accepts a signed byte cap. This helper validates
    the intended positive range ``1..127`` rather than reproducing Java's
    overflow and odd zero/negative-cap behavior for invalid parameters.
    """
    raw_mask = np.asarray(object_mask)
    if raw_mask.ndim != 3 or any(size == 0 for size in raw_mask.shape):
        raise ValueError("object_mask must be a nonempty 3D array in (x, y, z) order")
    if raw_mask.dtype != np.bool_:
        raise TypeError("object_mask must contain boolean values")

    try:
        cap = operator.index(max_distance)
    except TypeError as exc:
        raise TypeError("max_distance must be an integer") from exc
    if not 1 <= cap <= 127:
        raise ValueError("max_distance must be in the Java signed-byte range 1..127")

    mask = raw_mask
    distance = np.full(mask.shape, -1, dtype=np.int8)
    distance[mask] = 0

    frontier = np.zeros(mask.shape, dtype=bool)
    frontier[1:, :, :] |= mask[:-1, :, :]
    frontier[:-1, :, :] |= mask[1:, :, :]
    frontier[:, 1:, :] |= mask[:, :-1, :]
    frontier[:, :-1, :] |= mask[:, 1:, :]
    frontier[:, :, 1:] |= mask[:, :, :-1]
    frontier[:, :, :-1] |= mask[:, :, 1:]
    frontier &= ~mask

    distance[frontier] = 1

    # This is equivalent to the Java full-volume scans: each scan propagates
    # one distance layer. Java's first scan runs even for cap=1, producing the
    # source's distance-2 quirk at that cap.
    for layer in range(2, max(cap, 2) + 1):
        next_frontier = np.zeros(mask.shape, dtype=bool)
        next_frontier[1:, :, :] |= frontier[:-1, :, :]
        next_frontier[:-1, :, :] |= frontier[1:, :, :]
        next_frontier[:, 1:, :] |= frontier[:, :-1, :]
        next_frontier[:, :-1, :] |= frontier[:, 1:, :]
        next_frontier[:, :, 1:] |= frontier[:, :, :-1]
        next_frontier[:, :, :-1] |= frontier[:, :, 1:]
        next_frontier &= distance < 0
        if not np.any(next_frontier):
            break
        distance[next_frontier] = layer
        frontier = next_frontier

    return distance
