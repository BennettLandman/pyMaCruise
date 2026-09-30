"""TOADS relation-neighborhood maps from the selected VandyCRUISE source."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


# Source nested i/j/l loops visit these six faces in this order.
_FACES = ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0))
_INTERIOR = (slice(1, -1), slice(1, -1), slice(1, -1))
_EMPTY = np.int8(-1)


def toads_closest_labels(segmentation: ArrayLike, limit: int) -> NDArray[np.int8]:
    """Build the selected source's four-slot closest-label map.

    This implements ``MultipleToadSegmentation.computeClosestLabels``. Slot
    zero is the voxel's own class, slots one through three are distinct face
    neighbor labels, and missing entries are Java's ``EMPTY`` byte (-1).
    Missing slots diffuse synchronously for source passes ``1 <= t < 5*limit``.
    The selected method scans only interior voxels while seeding/diffusing;
    initialized border slots therefore remain unchanged. Operations are
    vectorized over voxels while preserving the Java face traversal order.
    """
    labels = np.asarray(segmentation)
    if labels.ndim != 3 or any(size == 0 for size in labels.shape):
        raise ValueError("segmentation must be a nonempty 3-D integer array")
    if not np.issubdtype(labels.dtype, np.integer) or np.issubdtype(labels.dtype, np.bool_):
        raise ValueError("segmentation labels must be integers")
    if any(size < 3 for size in labels.shape):
        raise ValueError("segmentation dimensions must be at least three voxels")
    if np.any(labels < 0) or np.any(labels > 127):
        raise ValueError("segmentation labels must fit the nonnegative Java byte class range")
    if isinstance(limit, (bool, np.bool_)) or not isinstance(limit, (int, np.integer)) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")

    class_labels = labels.astype(np.int8, copy=False)
    closest = np.full((*labels.shape, 4), _EMPTY, dtype=np.int8)
    closest[..., 0] = class_labels
    interior = closest[_INTERIOR]

    # Seed each interior voxel's distinct face labels with the selected loop order.
    for dx, dy, dz in _FACES:
        neighbor = class_labels[
            slice(1 + dx, labels.shape[0] - 1 + dx),
            slice(1 + dy, labels.shape[1] - 1 + dy),
            slice(1 + dz, labels.shape[2] - 1 + dz),
        ]
        slot1, slot2, slot3 = (interior[..., slot] for slot in (1, 2, 3))
        different = neighbor != interior[..., 0]
        fill1 = (slot1 == _EMPTY) & different
        fill2 = (~(slot1 == _EMPTY)) & (slot2 == _EMPTY) & (slot1 != neighbor) & different
        fill3 = (
            (~fill1)
            & (~fill2)
            & different
            & (slot3 == _EMPTY)
            & (slot1 != neighbor)
            & (slot2 != neighbor)
        )
        slot1[fill1] = neighbor[fill1]
        slot2[fill2] = neighbor[fill2]
        slot3[fill3] = neighbor[fill3]

    for _ in range(1, 5 * int(limit)):
        updated = closest.copy()
        target = updated[_INTERIOR]
        current = closest
        needs_slot1 = target[..., 1] == _EMPTY
        needs_slot2 = (target[..., 1] != _EMPTY) & (target[..., 2] == _EMPTY)
        needs_slot3 = (
            (target[..., 1] != _EMPTY)
            & (target[..., 2] != _EMPTY)
            & (target[..., 3] == _EMPTY)
        )
        for dx, dy, dz in _FACES:
            neighbor = current[
                slice(1 + dx, labels.shape[0] - 1 + dx),
                slice(1 + dy, labels.shape[1] - 1 + dy),
                slice(1 + dz, labels.shape[2] - 1 + dz),
                :,
            ]
            slot1 = target[..., 1]
            slot2 = target[..., 2]
            slot3 = target[..., 3]
            n1, n2, n3 = neighbor[..., 1], neighbor[..., 2], neighbor[..., 3]
            accept1 = needs_slot1 & (n1 != _EMPTY) & (target[..., 0] != n1)
            accept21 = (
                needs_slot2
                & (n1 != _EMPTY)
                & (target[..., 0] != n1)
                & (slot1 != n1)
            )
            accept22 = (
                needs_slot2
                & ~((n1 != _EMPTY) & (target[..., 0] != n1) & (slot1 != n1))
                & (n2 != _EMPTY)
                & (target[..., 0] != n2)
                & (slot1 != n2)
            )
            accept31 = (
                needs_slot3
                & (n1 != _EMPTY)
                & (target[..., 0] != n1)
                & (slot1 != n1)
                & (slot2 != n1)
            )
            accept32 = (
                needs_slot3
                & ~((n1 != _EMPTY) & (target[..., 0] != n1) & (slot1 != n1) & (slot2 != n1))
                & (n2 != _EMPTY)
                & (target[..., 0] != n2)
                & (slot1 != n2)
                & (slot2 != n2)
            )
            accept33 = (
                needs_slot3
                & ~(
                    ((n1 != _EMPTY) & (target[..., 0] != n1) & (slot1 != n1) & (slot2 != n1))
                    | ((n2 != _EMPTY) & (target[..., 0] != n2) & (slot1 != n2) & (slot2 != n2))
                )
                & (n3 != _EMPTY)
                & (target[..., 0] != n3)
                & (slot1 != n3)
                & (slot2 != n3)
            )
            slot1[accept1] = n1[accept1]
            slot2[accept21] = n1[accept21]
            slot2[accept22] = n2[accept22]
            slot3[accept31] = n1[accept31]
            slot3[accept32] = n2[accept32]
            slot3[accept33] = n3[accept33]
        closest[_INTERIOR] = target

    return closest
