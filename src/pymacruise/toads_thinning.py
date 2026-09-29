"""Selected VandyCRUISE TOADS thinning operations.

Behavior follows ``MultipleToadSegmentation.java`` (April 2006, author
Pierre-Louis Bazin), especially ``initGrowingLabels``,
``initThinningLabels``, and ``propagateThinning``. See synthetic
Java-reference reports for source hashes, matched slices, and known limits;
this module is not a complete TOADS port.

Paper context: heap-ordered object thinning/growing and topology-constrained
label changes implement components of the multiple-object fast-marching and
topology-preserving strategy (Bazin & Pham, 2007, 2008). Queue transitions
follow the selected Java source; the current validated slices do not establish
full-volume convergence or full TOADS behavior. See ``docs/papers/toads.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.toads_heap import ToadsHeapItem, ToadsPriorityHeap4D
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
    processing_label: int | None = None,
    empty_label: int = -1,
) -> ToadsThinningTransition:
    """Apply one source-ordered ``propagateThinning`` queue-item transition.

    This covers local target selection, relation-connected fallback by the
    largest neighbor shape prior, the established 18/6 multi-label topology
    gate, and the source volume/skeleton/segmentation updates. It assumes the
    caller already performed global queue initialization and popped an item
    for the center voxel. ``processing_label`` is the popped queue's class k;
    it defaults to the center segmentation label for direct local use. It does
    not implement queue order or whole-volume thinning.
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
    segmentation_center_label = int(segmentation[1, 1, 1])
    old_label = segmentation_center_label if processing_label is None else int(processing_label)
    skeleton_label = int(skeleton[1, 1, 1])
    if not 0 <= old_label < classes:
        raise ValueError("processing_label must be a valid class index")
    if skeleton_label != old_label:
        return ToadsThinningTransition(
            segmentation_center_label, skeleton_label, updated_volumes,
            available.copy(), None, None, False,
        )

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

    updated_segmentation = candidate if change_boundary else segmentation_center_label
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


@dataclass(frozen=True)
class ToadsGrowingQueueInitialization:
    """Availability and min-heap entries emitted by source initGrowingLabels."""

    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    queue: ToadsPriorityHeap4D


def toads_initialize_growing_queue(
    skeleton: ArrayLike,
    memberships: ArrayLike,
    ordering: ArrayLike,
    *,
    initial_available: ArrayLike | None = None,
    empty_label: int = -1,
    zero_threshold: float = 1e-10,
) -> ToadsGrowingQueueInitialization:
    """Initialize the selected VandyCRUISE initGrowingLabels queue.

    Interior empty-skeleton voxels begin available for every class; nonempty
    skeleton voxels begin unavailable. A candidate is then inserted for each
    available class touching that class on a six-face, using source priority
    ordering + maxMembershipRatio in a min-heap. Existing border availability
    is preserved, but it must be false so the source neighbor reads remain
    in bounds.
    """
    skel = np.asarray(skeleton)
    mems = np.asarray(memberships, dtype=np.float32)
    orderings = np.asarray(ordering, dtype=np.float32)
    if skel.ndim != 3 or any(length < 3 for length in skel.shape):
        raise ValueError("skeleton must be a 3-D grid with every dimension at least 3")
    if not np.issubdtype(skel.dtype, np.integer):
        raise ValueError("skeleton labels must be integers")
    if mems.ndim != 4 or mems.shape[:3] != skel.shape or mems.shape[-1] < 2:
        raise ValueError("memberships must be 4-D, match skeleton, and contain at least two classes")
    classes = mems.shape[-1]
    if not np.isfinite(mems).all():
        raise ValueError("memberships must be finite")
    if orderings.shape != skel.shape or not np.isfinite(orderings).all():
        raise ValueError("ordering must be finite and match the skeleton grid")
    if np.any((skel != empty_label) & ((skel < 0) | (skel >= classes))):
        raise ValueError("skeleton labels must be valid class indices or empty_label")
    zero = np.float32(zero_threshold)
    if not np.isfinite(zero) or zero < np.float32(0.0):
        raise ValueError("zero_threshold must be finite and nonnegative")

    if initial_available is None:
        available = np.zeros(skel.shape + (classes,), dtype=np.bool_)
    else:
        source_available = np.asarray(initial_available)
        if source_available.shape != skel.shape + (classes,) or source_available.dtype != np.bool_:
            raise ValueError("initial_available must be boolean and match skeleton and classes")
        available = source_available.copy()
    border = np.zeros(skel.shape, dtype=np.bool_)
    for axis in range(3):
        selectors = [slice(None)] * 3
        selectors[axis] = 0
        border[tuple(selectors)] = True
        selectors[axis] = -1
        border[tuple(selectors)] = True
    if np.any(available[border]):
        raise ValueError("initial border availability must be false for source-safe six-face reads")

    nx, ny, nz = skel.shape
    updated_skeleton = skel.astype(np.int16, copy=True)
    queue = ToadsPriorityHeap4D(nx * ny + ny * nz + nz * nx, "min")
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                is_empty = int(skel[x, y, z]) == empty_label
                available[x, y, z, :] = is_empty

    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for label in range(classes):
                    if not available[x, y, z, label]:
                        continue
                    boundary = any(
                        int(updated_skeleton[x + dx, y + dy, z + dz]) == label
                        for dx, dy, dz in _FACE_NEIGHBORS
                    )
                    if not boundary:
                        continue
                    available[x, y, z, label] = False
                    best = 0
                    maximum = np.float32(0.0)
                    for candidate in range(classes):
                        if mems[x, y, z, candidate] > maximum:
                            best = candidate
                            maximum = mems[x, y, z, candidate]
                    if label == best:
                        ratio = np.float32(1.0)
                    else:
                        ratio = np.float32(
                            np.float32(mems[x, y, z, best] + zero)
                            / np.float32(mems[x, y, z, label] + zero)
                        )
                    priority = np.float32(orderings[x, y, z] + ratio)
                    queue.add(priority, x, y, z, label)

    return ToadsGrowingQueueInitialization(updated_skeleton, available, queue)


@dataclass(frozen=True)
class ToadsGrowingQueueStep:
    """Image state after one selected source propagateGrowing event."""

    segmentation: NDArray[np.int16]
    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    ordering: NDArray[np.float32]
    class_volumes: NDArray[np.float32]
    popped_item: ToadsHeapItem
    topology_accepted: bool | None
    inserted_neighbors: tuple[ToadsHeapItem, ...]


@dataclass(frozen=True)
class ToadsGrowingPass:
    """Final state and event trace from one source growing queue pass."""

    segmentation: NDArray[np.int16]
    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    ordering: NDArray[np.float32]
    class_volumes: NDArray[np.float32]
    events: tuple[ToadsGrowingQueueStep, ...]
    remaining_queue: tuple[ToadsHeapItem, ...]


def toads_process_next_growing_item(
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    memberships: ArrayLike,
    available: ArrayLike,
    ordering: ArrayLike,
    class_volumes: ArrayLike,
    queue: ToadsPriorityHeap4D,
    *,
    empty_label: int = -1,
    zero_threshold: float = 1e-10,
) -> ToadsGrowingQueueStep:
    """Process one source-ordered topology-gated growing item."""
    labels = np.asarray(segmentation)
    skel = np.asarray(skeleton)
    mems = np.asarray(memberships, dtype=np.float32)
    avail = np.asarray(available)
    orderings = np.asarray(ordering, dtype=np.float32)
    volumes = np.asarray(class_volumes, dtype=np.float32)
    if labels.ndim != 3 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be a 3-D integer grid")
    if skel.shape != labels.shape or not np.issubdtype(skel.dtype, np.integer):
        raise ValueError("skeleton must be an integer grid matching segmentation")
    if mems.ndim != 4 or mems.shape[:3] != labels.shape or mems.shape[-1] < 2:
        raise ValueError("memberships must be 4-D and match segmentation")
    classes = mems.shape[-1]
    if avail.shape != labels.shape + (classes,) or avail.dtype != np.bool_:
        raise ValueError("available must be a boolean grid matching segmentation and classes")
    if orderings.shape != labels.shape or volumes.shape != (classes,):
        raise ValueError("ordering and class volumes must match the image and class count")
    if not np.isfinite(mems).all() or not np.isfinite(orderings).all() or not np.isfinite(volumes).all():
        raise ValueError("memberships, ordering, and class volumes must be finite")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation labels must be valid class indices")
    if np.any((skel != empty_label) & ((skel < 0) | (skel >= classes))):
        raise ValueError("skeleton labels must be valid class indices or empty_label")
    zero = np.float32(zero_threshold)
    if not np.isfinite(zero) or zero < np.float32(0.0):
        raise ValueError("zero_threshold must be finite and nonnegative")
    if not isinstance(queue, ToadsPriorityHeap4D) or queue.mode != "min":
        raise ValueError("queue must be a TOADS min heap")
    if not queue:
        raise ValueError("cannot process an empty TOADS growing queue")

    next_item = queue.peek_first()
    x, y, z, label = next_item.x, next_item.y, next_item.z, next_item.label
    if not (
        1 <= x < labels.shape[0] - 1
        and 1 <= y < labels.shape[1] - 1
        and 1 <= z < labels.shape[2] - 1
        and 0 <= label < classes
    ):
        raise ValueError("growing queue item must have an interior coordinate and valid class")
    popped = queue.pop_first()
    updated_labels = labels.astype(np.int16, copy=True)
    updated_skeleton = skel.astype(np.int16, copy=True)
    updated_available = avail.copy()
    updated_ordering = orderings.copy()
    updated_volumes = volumes.copy()
    topology_accepted: bool | None = None
    inserted: list[ToadsHeapItem] = []

    if int(updated_skeleton[x, y, z]) == empty_label:
        neighborhood = updated_labels[x - 1 : x + 2, y - 1 : y + 2, z - 1 : z + 2]
        topology_accepted = toads_homeomorphic_label_change(neighborhood, label)
        if topology_accepted:
            updated_skeleton[x, y, z] = np.int16(label)
            updated_available[x, y, z, :] = False
            updated_ordering[x, y, z] = popped.value
            updated_volumes[label] = np.float32(volumes[label] + mems[x, y, z, label])
            updated_labels[x, y, z] = np.int16(label)
            for dx, dy, dz in _FACE_NEIGHBORS:
                xi, yj, zl = x + dx, y + dy, z + dz
                if not updated_available[xi, yj, zl, label]:
                    continue
                best = 0
                maximum = np.float32(0.0)
                for candidate in range(classes):
                    if mems[xi, yj, zl, candidate] > maximum:
                        best = candidate
                        maximum = mems[xi, yj, zl, candidate]
                if label == best:
                    ratio = np.float32(1.0)
                else:
                    ratio = np.float32(
                        np.float32(mems[xi, yj, zl, best] + zero)
                        / np.float32(mems[xi, yj, zl, label] + zero)
                    )
                priority = np.float32(popped.value + ratio)
                updated_available[xi, yj, zl, label] = False
                queue.add(priority, xi, yj, zl, label)
                inserted.append(ToadsHeapItem(priority, xi, yj, zl, label))
        else:
            updated_available[x, y, z, label] = True

    return ToadsGrowingQueueStep(
        updated_labels,
        updated_skeleton,
        updated_available,
        updated_ordering,
        updated_volumes,
        popped,
        topology_accepted,
        tuple(inserted),
    )


def toads_propagate_growing(
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    memberships: ArrayLike,
    ordering: ArrayLike,
    class_volumes: ArrayLike,
    *,
    initial_available: ArrayLike | None = None,
    empty_label: int = -1,
    zero_threshold: float = 1e-10,
) -> ToadsGrowingPass:
    """Initialize and exhaust the selected ``propagateGrowing`` min-heap.

    This composes ``initGrowingLabels`` and one-event growing transitions until
    the source ``while (isNotEmpty())`` loop terminates. It returns the full
    event trace for deterministic source comparison; it is not full TOADS
    segmentation or atlas registration.
    """
    labels = np.asarray(segmentation)
    source_skeleton = np.asarray(skeleton)
    mems = np.asarray(memberships, dtype=np.float32)
    orderings = np.asarray(ordering, dtype=np.float32)
    volumes = np.asarray(class_volumes, dtype=np.float32)
    if labels.ndim != 3 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be a 3-D integer grid")
    if source_skeleton.shape != labels.shape:
        raise ValueError("skeleton must match segmentation")
    if mems.ndim != 4 or mems.shape[:3] != labels.shape:
        raise ValueError("memberships must be 4-D and match segmentation")
    if np.any(labels < 0) or np.any(labels >= mems.shape[-1]):
        raise ValueError("segmentation labels must be valid membership class indices")
    if volumes.shape != (mems.shape[-1],) or not np.isfinite(volumes).all():
        raise ValueError("class_volumes must be finite and match membership classes")
    initialized = toads_initialize_growing_queue(
        source_skeleton,
        memberships,
        ordering,
        initial_available=initial_available,
        empty_label=empty_label,
        zero_threshold=zero_threshold,
    )
    current_segmentation = labels.astype(np.int16, copy=True)
    current_skeleton = initialized.skeleton
    current_available = initialized.available
    current_ordering = orderings.copy()
    current_volumes = volumes.copy()
    events: list[ToadsGrowingQueueStep] = []
    while initialized.queue:
        event = toads_process_next_growing_item(
            current_segmentation,
            current_skeleton,
            memberships,
            current_available,
            current_ordering,
            current_volumes,
            initialized.queue,
            empty_label=empty_label,
            zero_threshold=zero_threshold,
        )
        events.append(event)
        current_segmentation = event.segmentation
        current_skeleton = event.skeleton
        current_available = event.available
        current_ordering = event.ordering
        current_volumes = event.class_volumes

    remaining: list[ToadsHeapItem] = []
    while initialized.queue:
        remaining.append(initialized.queue.pop_first())
    return ToadsGrowingPass(
        current_segmentation,
        current_skeleton,
        current_available,
        current_ordering,
        current_volumes,
        tuple(events),
        tuple(remaining),
    )


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


def toads_schedule_thinning_neighbors(
    segmentation: ArrayLike,
    memberships: ArrayLike,
    available: ArrayLike,
    queue: ToadsPriorityHeap4D,
    removed_voxel: tuple[int, int, int],
    removed_label: int,
    popped_priority: float,
    class_volume: float,
    minimum_class_volume: float,
    limit: float,
    *,
    zero_threshold: float = 1e-10,
) -> tuple[NDArray[np.bool_], tuple[ToadsHeapItem, ...]]:
    """Enqueue source-eligible face neighbors after a thinning removal.

    This is the post-removal neighbor loop in ``propagateThinning``. Callers
    supply the already-updated segmentation, heap, and class volume. The
    source code's inner neighbor count tests the *outer* face-offset norm,
    which is always nonzero here; that literal behavior is intentionally
    retained for reference matching.
    """
    labels = np.asarray(segmentation)
    mems = np.asarray(memberships, dtype=np.float32)
    state = np.asarray(available)
    if labels.ndim != 3 or any(length < 3 for length in labels.shape):
        raise ValueError("segmentation must be a 3-D grid with every dimension at least 3")
    if not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation labels must be integers")
    if mems.ndim != 4 or mems.shape[:3] != labels.shape or mems.shape[-1] < 2:
        raise ValueError("memberships must be 4-D and match segmentation")
    classes = mems.shape[-1]
    if state.shape != labels.shape + (classes,) or state.dtype != np.bool_:
        raise ValueError("available must be a boolean grid matching segmentation and classes")
    if not np.isfinite(mems).all():
        raise ValueError("memberships must be finite")
    if not isinstance(queue, ToadsPriorityHeap4D) or queue.mode != "max":
        raise ValueError("queue must be a TOADS max heap")
    if not isinstance(removed_label, (int, np.integer)) or not 0 <= removed_label < classes:
        raise ValueError("removed_label must be a valid class index")
    x, y, z = removed_voxel
    if not (1 <= x < labels.shape[0] - 1 and 1 <= y < labels.shape[1] - 1 and 1 <= z < labels.shape[2] - 1):
        raise ValueError("removed_voxel must be interior")
    zero = np.float32(zero_threshold)
    value = np.float32(popped_priority)
    volume = np.float32(class_volume)
    minimum = np.float32(minimum_class_volume)
    limit32 = np.float32(limit)
    if not all(np.isfinite(v) for v in (zero, value, volume, minimum, limit32)) or zero < 0.0 or limit32 < 0.0:
        raise ValueError("priority, volume, limit, and ZERO must be finite; ZERO/limit nonnegative")

    updated_available = state.copy()
    inserted = []
    # Same nested i/j/l order as BinaryHeap4D's selected caller.
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                if i * i + j * j + k * k != 1:
                    continue
                xi, yj, zl = x + i, y + j, z + k
                if not updated_available[xi, yj, zl, removed_label]:
                    continue
                neighbor_count = 0
                for ni in (-1, 0, 1):
                    for nj in (-1, 0, 1):
                        for nl in (-1, 0, 1):
                            # Source quirk: this condition refers to outer
                            # i/j/k, not inner ni/nj/nl, so it is always true.
                            if i * i + j * j + k * k > 0 and labels[xi + ni, yj + nj, zl + nl] == removed_label:
                                neighbor_count += 1
                if neighbor_count == 1:
                    order = value
                else:
                    best = 0
                    maximum = np.float32(0.0)
                    for candidate in range(classes):
                        if mems[xi, yj, zl, candidate] > maximum:
                            best = candidate
                            maximum = mems[xi, yj, zl, candidate]
                    if removed_label == best:
                        ratio = np.float32(1.0)
                    else:
                        ratio = np.float32(
                            np.float32(mems[xi, yj, zl, best] + zero)
                            / np.float32(mems[xi, yj, zl, removed_label] + zero)
                        )
                    order = np.float32(value - np.float32(np.float32(1.0) / ratio))
                if volume > minimum and order > np.float32(-limit32):
                    updated_available[xi, yj, zl, removed_label] = False
                    queue.add(order, xi, yj, zl, removed_label)
                    inserted.append(ToadsHeapItem(np.float32(order), xi, yj, zl, removed_label))

    return updated_available, tuple(inserted)


@dataclass(frozen=True)
class ToadsThinningQueueStep:
    """Full arrays and local decision after processing one max-heap item."""

    segmentation: NDArray[np.int16]
    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    class_volumes: NDArray[np.float32]
    popped_item: ToadsHeapItem
    transition: ToadsThinningTransition | None
    inserted_neighbors: tuple[ToadsHeapItem, ...]


@dataclass(frozen=True)
class ToadsThinningQueueRun:
    """Final state after the source outer loop drains or crosses its limit."""

    segmentation: NDArray[np.int16]
    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    class_volumes: NDArray[np.float32]
    processed_items: tuple[ToadsHeapItem, ...]
    queue: ToadsPriorityHeap4D


@dataclass(frozen=True)
class ToadsThinningPass:
    """Initialized source thinning state after one complete queue pass."""

    segmentation: NDArray[np.int16]
    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    class_volumes: NDArray[np.float32]
    boundary_counts: NDArray[np.int32]
    minimum_class_volumes: NDArray[np.float32]
    processed_items: tuple[ToadsHeapItem, ...]
    queue: ToadsPriorityHeap4D


def toads_propagate_thinning(
    segmentation: ArrayLike,
    memberships: ArrayLike,
    relation_matrix: ArrayLike,
    transformed_shape_values: ArrayLike,
    volume_ratio: float,
    limit: float,
    *,
    empty_label: int = -1,
    zero_threshold: float = 1e-10,
) -> ToadsThinningPass:
    """Run source ``initThinningLabels`` and one complete thinning pass.

    The source fixes its absolute minimum class volume at 100 and combines it
    with ``volume_ratio * initial_volume`` using Java float32 arithmetic.
    This composes only the exact-membership thinning path; growing and the
    larger TOADS iteration remain separate.
    """
    ratio32 = np.float32(volume_ratio)
    if not np.isfinite(ratio32):
        raise ValueError("volume_ratio must be finite and representable as float32")
    initialized = toads_initialize_thinning_queue(
        segmentation,
        memberships,
        empty_label=empty_label,
        zero_threshold=zero_threshold,
    )
    minimums = np.maximum(
        np.float32(100.0), np.float32(ratio32 * initialized.class_volumes)
    ).astype(np.float32, copy=False)
    run = toads_run_thinning_queue(
        segmentation,
        initialized.skeleton,
        memberships,
        initialized.available,
        relation_matrix,
        transformed_shape_values,
        initialized.class_volumes,
        minimums,
        initialized.queue,
        limit,
        empty_label=empty_label,
        zero_threshold=zero_threshold,
    )
    return ToadsThinningPass(
        run.segmentation,
        run.skeleton,
        run.available,
        run.class_volumes,
        initialized.boundary_counts,
        minimums,
        run.processed_items,
        run.queue,
    )


def toads_run_thinning_queue(
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    memberships: ArrayLike,
    available: ArrayLike,
    relation_matrix: ArrayLike,
    transformed_shape_values: ArrayLike,
    class_volumes: ArrayLike,
    minimum_class_volumes: ArrayLike,
    queue: ToadsPriorityHeap4D,
    limit: float,
    *,
    empty_label: int = -1,
    zero_threshold: float = 1e-10,
) -> ToadsThinningQueueRun:
    """Run the selected ``propagateThinning`` outer heap loop.

    The source initializes ``val`` to 0 and tests its previous popped value in
    ``while (isNotEmpty() && (val > -limit))``. Consequently, the first item
    at or below the negative limit is still popped and processed if the prior
    item was above it. This composes the queue-event port only; it does not
    initialize thinning labels or perform the rest of TOADS.
    """
    limit32 = np.float32(limit)
    if not np.isfinite(limit32) or limit32 < np.float32(0.0):
        raise ValueError("limit must be finite, nonnegative, and representable as float32")
    if not isinstance(queue, ToadsPriorityHeap4D) or queue.mode != "max":
        raise ValueError("queue must be a TOADS max heap")
    labels = np.asarray(segmentation, dtype=np.int16).copy()
    skel = np.asarray(skeleton, dtype=np.int16).copy()
    available_state = np.asarray(available, dtype=np.bool_).copy()
    volumes = np.asarray(class_volumes, dtype=np.float32).copy()
    previous_priority = np.float32(0.0)
    processed: list[ToadsHeapItem] = []
    while queue and previous_priority > np.float32(-limit32):
        step = toads_process_next_thinning_item(
            labels,
            skel,
            memberships,
            available_state,
            relation_matrix,
            transformed_shape_values,
            volumes,
            minimum_class_volumes,
            queue,
            limit32,
            empty_label=empty_label,
            zero_threshold=zero_threshold,
        )
        labels = step.segmentation
        skel = step.skeleton
        available_state = step.available
        volumes = step.class_volumes
        previous_priority = step.popped_item.value
        processed.append(step.popped_item)
    return ToadsThinningQueueRun(
        labels, skel, available_state, volumes, tuple(processed), queue
    )


def toads_process_next_thinning_item(
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    memberships: ArrayLike,
    available: ArrayLike,
    relation_matrix: ArrayLike,
    transformed_shape_values: ArrayLike,
    class_volumes: ArrayLike,
    minimum_class_volumes: ArrayLike,
    queue: ToadsPriorityHeap4D,
    limit: float,
    *,
    empty_label: int = -1,
    zero_threshold: float = 1e-10,
) -> ToadsThinningQueueStep:
    """Process one selected `propagateThinning` heap item and its neighbors.

    This composes the independently validated queue, local topology/volume
    transition, and post-removal six-face scheduling pieces. Inputs are copied;
    the supplied queue is popped and receives any newly scheduled neighbors.
    A stale item whose label no longer matches the current skeleton is
    discarded without changing image state. As in the source, the queue label
    selects the class volume and membership even if segmentation differs.
    """
    labels = np.asarray(segmentation)
    skel = np.asarray(skeleton)
    mems = np.asarray(memberships, dtype=np.float32)
    avail = np.asarray(available)
    relations = np.asarray(relation_matrix)
    shapes = np.asarray(transformed_shape_values, dtype=np.float32)
    volumes = np.asarray(class_volumes, dtype=np.float32)
    minimums = np.asarray(minimum_class_volumes, dtype=np.float32)
    if labels.ndim != 3 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be a 3-D integer grid")
    if skel.shape != labels.shape or not np.issubdtype(skel.dtype, np.integer):
        raise ValueError("skeleton must be an integer grid matching segmentation")
    if mems.ndim != 4 or mems.shape[:3] != labels.shape or mems.shape[-1] < 2:
        raise ValueError("memberships must be 4-D and match segmentation")
    classes = mems.shape[-1]
    if avail.shape != labels.shape + (classes,) or avail.dtype != np.bool_:
        raise ValueError("available must be a boolean grid matching segmentation and classes")
    if relations.shape != (classes, classes) or relations.dtype != np.bool_:
        raise ValueError("relation_matrix must be a boolean square matrix matching classes")
    if shapes.shape != mems.shape or not np.isfinite(shapes).all():
        raise ValueError("transformed_shape_values must be finite and match memberships")
    if volumes.shape != (classes,) or minimums.shape != (classes,):
        raise ValueError("class volume arrays must match class count")
    if not np.isfinite(volumes).all() or not np.isfinite(minimums).all() or not np.isfinite(mems).all():
        raise ValueError("memberships and class volumes must be finite")
    if not isinstance(queue, ToadsPriorityHeap4D) or queue.mode != "max":
        raise ValueError("queue must be a TOADS max heap")
    if not queue:
        raise ValueError("cannot process an empty TOADS thinning queue")

    updated_labels = labels.astype(np.int16, copy=True)
    updated_skeleton = skel.astype(np.int16, copy=True)
    updated_available = avail.copy()
    updated_volumes = volumes.copy()
    popped = queue.pop_first()
    x, y, z, label = popped.x, popped.y, popped.z, popped.label
    transition: ToadsThinningTransition | None = None
    inserted: tuple[ToadsHeapItem, ...] = ()
    if (
        1 <= x < labels.shape[0] - 1
        and 1 <= y < labels.shape[1] - 1
        and 1 <= z < labels.shape[2] - 1
        and 0 <= label < classes
        and int(updated_skeleton[x, y, z]) == label
        and updated_volumes[label] > minimums[label]
    ):
        transition = toads_thinning_voxel_transition(
            updated_labels[x - 1 : x + 2, y - 1 : y + 2, z - 1 : z + 2],
            updated_skeleton[x - 1 : x + 2, y - 1 : y + 2, z - 1 : z + 2],
            mems[x, y, z, :],
            relations,
            shapes[x - 1 : x + 2, y - 1 : y + 2, z - 1 : z + 2, :],
            updated_volumes,
            minimums,
            available_at_center=updated_available[x, y, z, :],
            processing_label=label,
            empty_label=empty_label,
        )
        updated_labels[x, y, z] = transition.segmentation_label
        updated_skeleton[x, y, z] = transition.skeleton_label
        updated_available[x, y, z, :] = transition.available_at_center
        updated_volumes = transition.class_volumes.copy()
        if transition.eligible:
            updated_available, inserted = toads_schedule_thinning_neighbors(
                updated_labels,
                mems,
                updated_available,
                queue,
                (x, y, z),
                label,
                popped.value,
                updated_volumes[label],
                minimums[label],
                limit,
                zero_threshold=zero_threshold,
            )
    return ToadsThinningQueueStep(
        updated_labels, updated_skeleton, updated_available, updated_volumes,
        popped, transition, inserted,
    )
