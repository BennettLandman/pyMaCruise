"""TOADS relation-neighborhood maps from the selected VandyCRUISE source."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


# Source nested i/j/l loops visit these six faces in this order.
_FACES = ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0))
_INTERIOR = (slice(1, -1), slice(1, -1), slice(1, -1))
_EMPTY = np.int8(-1)


@dataclass(frozen=True)
class ToadsRelationList:
    """Source-style relation map and relation-membership table."""

    relation_map: NDArray[np.int16]
    relation_membership: NDArray[np.bool_]

    @property
    def class_adjacency(self) -> NDArray[np.bool_]:
        """Return the selected source's first ``class_count`` union rows."""
        return self.relation_membership[: self.relation_membership.shape[1]]

    @property
    def class_count(self) -> int:
        return int(self.relation_membership.shape[1])


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


def toads_build_relation_list(
    segmentation: ArrayLike,
    class_count: int,
) -> ToadsRelationList:
    """Build source-ordered relation IDs for local six-face label sets.

    This ports ``MultipleToadSegmentation.buildRelationList`` from a fresh
    constructor-equivalent relation map. Singleton relations occupy IDs
    ``0..class_count-1``. Each distinct boundary set, containing the voxel's
    own class and its unique six-face neighbors, receives the next ID in the
    source x/y/z first-occurrence order; set equality is order-independent.
    The returned table contains those sets followed by the source's unioned
    adjacency rows for each singleton class. Borders retain their initialized
    class labels, as the Java method writes only interior voxels. A uint64
    bitset is used internally to vectorize boundary set collection.
    """
    labels = np.asarray(segmentation)
    if labels.ndim != 3 or any(size == 0 for size in labels.shape):
        raise ValueError("segmentation must be a nonempty 3-D integer array")
    if not np.issubdtype(labels.dtype, np.integer) or np.issubdtype(labels.dtype, np.bool_):
        raise ValueError("segmentation labels must be integers")
    if isinstance(class_count, (bool, np.bool_)) or not isinstance(class_count, (int, np.integer)):
        raise ValueError("class_count must be a positive integer")
    classes = int(class_count)
    if not 1 <= classes <= 63:
        raise ValueError("class_count must be between 1 and 63 for source class bitsets")
    if any(size < 3 for size in labels.shape):
        raise ValueError("segmentation dimensions must be at least three voxels")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation contains an out-of-range class index")

    class_labels = labels.astype(np.int64, copy=False)
    center = class_labels[_INTERIOR]
    relation_keys = np.left_shift(np.uint64(1), center.astype(np.uint64))
    is_boundary = np.zeros(center.shape, dtype=bool)
    for dx, dy, dz in _FACES:
        neighbor = class_labels[
            slice(1 + dx, labels.shape[0] - 1 + dx),
            slice(1 + dy, labels.shape[1] - 1 + dy),
            slice(1 + dz, labels.shape[2] - 1 + dz),
        ]
        differs = neighbor != center
        is_boundary |= differs
        bits = np.left_shift(np.uint64(1), neighbor.astype(np.uint64))
        relation_keys |= np.where(differs, bits, np.uint64(0))

    boundary_keys = relation_keys[is_boundary]
    if boundary_keys.size:
        unique_keys, first_indices = np.unique(boundary_keys, return_index=True)
        ordered_keys = unique_keys[np.argsort(first_indices)]
    else:
        ordered_keys = np.empty(0, dtype=np.uint64)
    if classes + ordered_keys.size > 500:
        raise ValueError("source relation table exceeds the 500-entry capacity")

    relation_map = class_labels.astype(np.int16, copy=True)
    relation_map_interior = relation_map[_INTERIOR]
    for relation_id, key in enumerate(ordered_keys, start=classes):
        relation_map_interior[is_boundary & (relation_keys == key)] = relation_id

    relation_membership = np.zeros((classes + ordered_keys.size, classes), dtype=bool)
    relation_membership[np.arange(classes), np.arange(classes)] = True
    for index, key in enumerate(ordered_keys, start=classes):
        members = np.array(
            [bool((int(key) >> class_index) & 1) for class_index in range(classes)],
            dtype=bool,
        )
        relation_membership[index] = members
        for class_index in np.flatnonzero(members):
            relation_membership[class_index] |= members

    return ToadsRelationList(relation_map, relation_membership)
