"""Local selected VandyCRUISE TOADS thinning state transition."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.toads_heap import ToadsPriorityHeap4D
from pymacruise.toads_topology import toads_homeomorphic_label_change


_FACE_NEIGHBORS = ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0))


@dataclass(frozen=True)
class ToadsThinningTransition:
    """Discrete and volume outputs for one eligible thinning queue item."""

    segmentation_label: int
    skeleton_label: int
    class_volumes: NDArray[np.float32]
    available_at_center: NDArray[np.bool_]
    candidate_label: int | None
    topology_accepted: bool | None
    eligible: bool


def toads_thinning_voxel_transition(
    segmentation_neighborhood: ArrayLike,
    skeleton_neighborhood: ArrayLike,
    memberships: ArrayLike,
    relation_matrix: ArrayLike,
    transformed_shape_neighborhood: ArrayLike,
    class_volumes: ArrayLike,
    minimum_class_volumes: ArrayLike,
    *,
    available_at_center: ArrayLike | None = None,
    empty_label: int = -1,
) -> ToadsThinningTransition:
    """Apply one source-ordered ``propagateThinning`` queue-item transition.

    This covers local target selection, relation-connected fallback by the
    largest neighbor shape prior, the established 18/6 multi-label topology
    gate, and the source volume/skeleton/segmentation updates. It assumes the
    caller already performed global queue initialization and popped an item
    for the center voxel. It does not implement queue order or whole-volume
    thinning.
    """
    segmentation = np.asarray(segmentation_neighborhood)
    skeleton = np.asarray(skeleton_neighborhood)
    membership = np.asarray(memberships, dtype=np.float32)
    relations = np.asarray(relation_matrix)
    shapes = np.asarray(transformed_shape_neighborhood, dtype=np.float32)
    volumes = np.asarray(class_volumes, dtype=np.float32)
    minimums = np.asarray(minimum_class_volumes, dtype=np.float32)
    available = (np.zeros(membership.shape, dtype=np.bool_) if available_at_center is None
                 else np.asarray(available_at_center, dtype=np.bool_))
    if segmentation.shape != (3, 3, 3) or skeleton.shape != (3, 3, 3):
        raise ValueError("segmentation and skeleton neighborhoods must have shape (3, 3, 3)")
    if not np.issubdtype(segmentation.dtype, np.integer) or not np.issubdtype(skeleton.dtype, np.integer):
        raise ValueError("segmentation and skeleton labels must be integers")
    if not isinstance(empty_label, (int, np.integer)):
        raise ValueError("empty_label must be an integer")
    classes = membership.size
    if membership.ndim != 1 or classes < 2 or not np.isfinite(membership).all():
        raise ValueError("memberships must be a finite class vector with at least two classes")
    if relations.shape != (classes, classes) or relations.dtype != np.bool_:
        raise ValueError("relation_matrix must be a boolean square matrix matching classes")
    if shapes.shape != (3, 3, 3, classes) or not np.isfinite(shapes).all():
        raise ValueError("transformed_shape_neighborhood must be finite with shape (3, 3, 3, classes)")
    if volumes.shape != (classes,) or minimums.shape != (classes,):
        raise ValueError("class volume arrays must match the membership class count")
    if available.shape != (classes,):
        raise ValueError("available_at_center must match the membership class count")
    if not np.isfinite(volumes).all() or not np.isfinite(minimums).all():
        raise ValueError("class volumes and minimums must be finite")
    if np.any(segmentation < 0) or np.any(segmentation >= classes):
        raise ValueError("segmentation labels must be valid class indices")
    if np.any((skeleton != empty_label) & ((skeleton < 0) | (skeleton >= classes))):
        raise ValueError("skeleton labels must be valid class indices or empty_label")

    updated_volumes = volumes.copy()
    old_label = int(segmentation[1, 1, 1])
    skeleton_label = int(skeleton[1, 1, 1])
    if skeleton_label != old_label:
        return ToadsThinningTransition(old_label, skeleton_label, updated_volumes, available.copy(), None, None, False)

    # Source bestMembership uses strict greater-than, so ties stay with the
    # first class and all-nonpositive scores select class 0.
    target = 0
    maximum = np.float32(0.0)
    for label in range(classes):
        if membership[label] > maximum:
            target = label
            maximum = membership[label]

    candidate = target
    change_boundary = False
    if target != old_label:
        direct = False
        next_label: int | None = None
        best_prior = np.float32(-0.1)
        for dx, dy, dz in _FACE_NEIGHBORS:
            neighbor_label = int(segmentation[1 + dx, 1 + dy, 1 + dz])
            if neighbor_label == target:
                direct = True
            elif neighbor_label != old_label and relations[neighbor_label, target]:
                shape_prior = shapes[1 + dx, 1 + dy, 1 + dz, target]
                if shape_prior > best_prior:
                    next_label = neighbor_label
                    best_prior = shape_prior
        if direct:
            change_boundary = True
        elif next_label is not None:
            candidate = next_label
            change_boundary = True

    topology_accepted: bool | None = None
    if change_boundary:
        topology_accepted = toads_homeomorphic_label_change(segmentation, candidate)
        change_boundary = topology_accepted

    eligible = bool(volumes[old_label] > minimums[old_label])
    proposed_label = candidate if target != old_label and (topology_accepted is not None) else None
    if not eligible:
        return ToadsThinningTransition(old_label, skeleton_label, updated_volumes, available.copy(), proposed_label, topology_accepted, False)

    updated_segmentation = candidate if change_boundary else old_label
    updated_volumes[old_label] = np.float32(volumes[old_label] - membership[old_label])
    updated_available = available.copy()
    updated_available[old_label] = False
    return ToadsThinningTransition(
        updated_segmentation,
        int(empty_label),
        updated_volumes,
        updated_available,
        proposed_label,
        topology_accepted,
        True,
    )


@dataclass(frozen=True)
class ToadsThinningQueueInitialization:
    """Arrays and priority queue emitted by one source thinning initializer."""

    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    class_volumes: NDArray[np.float32]
    boundary_counts: NDArray[np.int32]
    queue: ToadsPriorityHeap4D


def toads_initialize_thinning_queue(
    segmentation: ArrayLike,
    memberships: ArrayLike,
    *,
    empty_label: int = -1,
    initial_skeleton: ArrayLike | None = None,
    initial_available: ArrayLike | None = None,
    zero_threshold: float = 1e-10,
) -> ToadsThinningQueueInitialization:
    """Initialize selected ``initThinningLabels`` state and max-heap entries.

    The source scans only interior voxels, accumulates current-class
    memberships in x/y/z order, detects boundaries from six face neighbors,
    retains the source endpoint test's neighbor-count behavior, and inserts
    negative membership-ratio priorities into a max heap. Existing border
    skeleton/availability values are copied through unchanged.
    """
    labels = np.asarray(segmentation)
    mems = np.asarray(memberships, dtype=np.float32)
    if labels.ndim != 3 or any(length < 3 for length in labels.shape):
        raise ValueError("segmentation must be a 3-D grid with every dimension at least 3")
    if not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation labels must be integers")
    if mems.ndim != 4 or mems.shape[:3] != labels.shape or mems.shape[-1] < 2:
        raise ValueError("memberships must be 4-D, match segmentation, and contain at least two classes")
    classes = mems.shape[-1]
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation labels must be valid class indices")
    if not np.isfinite(mems).all():
        raise ValueError("memberships must be finite")
    zero = np.float32(zero_threshold)
    if not np.isfinite(zero) or zero < np.float32(0.0):
        raise ValueError("zero_threshold must be finite and nonnegative")

    if initial_skeleton is None:
        skeleton = np.full(labels.shape, empty_label, dtype=np.int16)
    else:
        source_skeleton = np.asarray(initial_skeleton)
        if source_skeleton.shape != labels.shape or not np.issubdtype(source_skeleton.dtype, np.integer):
            raise ValueError("initial_skeleton must be an integer grid matching segmentation")
        skeleton = source_skeleton.astype(np.int16, copy=True)
    if initial_available is None:
        available = np.zeros(labels.shape + (classes,), dtype=np.bool_)
    else:
        source_available = np.asarray(initial_available)
        if source_available.shape != labels.shape + (classes,) or source_available.dtype != np.bool_:
            raise ValueError("initial_available must be a boolean grid matching segmentation and classes")
        available = source_available.copy()

    nx, ny, nz = labels.shape
    volumes = np.zeros(classes, dtype=np.float32)
    boundaries = np.zeros(classes, dtype=np.int32)
    queue = ToadsPriorityHeap4D(nx * ny + ny * nz + nz * nx, "max")

    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                current = int(labels[x, y, z])
                skeleton[x, y, z] = np.int16(current)
                volumes[current] = np.float32(volumes[current] + mems[x, y, z, current])
                for candidate in range(classes):
                    available[x, y, z, candidate] = candidate == current

    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                current = int(labels[x, y, z])
                for candidate in range(classes):
                    if current != candidate:
                        continue
                    boundary = False
                    neighbor_count = 0
                    ni = nj = nl = 0
                    for i in (-1, 0, 1):
                        for j in (-1, 0, 1):
                            for k in (-1, 0, 1):
                                distance2 = i * i + j * j + k * k
                                neighbor_label = int(labels[x + i, y + j, z + k])
                                if distance2 == 1 and neighbor_label != candidate:
                                    boundary = True
                                elif distance2 > 0:
                                    neighbor_count += 1
                                    ni, nj, nl = x + i, y + j, z + k
                    if not boundary:
                        continue
                    available[x, y, z, candidate] = False

                    endpoint = False
                    if neighbor_count == 1:
                        neighbor_count = 0
                        for i in (-1, 0, 1):
                            for j in (-1, 0, 1):
                                for k in (-1, 0, 1):
                                    if i * i + j * j + k * k > 0 and labels[ni + i, nj + j, nl + k] == candidate:
                                        neighbor_count += 1
                        endpoint = neighbor_count <= 2

                    if endpoint:
                        order = np.float32(-0.0)
                    else:
                        best = 0
                        maximum = np.float32(0.0)
                        for other in range(classes):
                            if mems[x, y, z, other] > maximum:
                                best = other
                                maximum = mems[x, y, z, other]
                        if candidate == best:
                            ratio = np.float32(1.0)
                        else:
                            ratio = np.float32(
                                np.float32(mems[x, y, z, best] + zero)
                                / np.float32(mems[x, y, z, candidate] + zero)
                            )
                        order = np.float32(np.float32(-1.0) / ratio)
                    queue.add(order, x, y, z, candidate)
                    boundaries[candidate] += 1

    return ToadsThinningQueueInitialization(skeleton, available, volumes, boundaries, queue)
