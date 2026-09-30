"""Selected LesionTOADS growing-frontier operations.

This ports the selected VandyCRUISE ``LesionToadSegmentation`` methods
``initGrowingLabels``, ``bestMembership``, and ``maxMembershipRatio`` through
queue exhaustion. It reuses the already source-validated ``BinaryHeap4D``
port in :mod:`pymacruise.toads_heap`; it does not estimate lesion memberships
or run later competition/thinning stages.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.toads_heap import ToadsHeapItem, ToadsPriorityHeap4D
from pymacruise.toads_topology import toads_homeomorphic_label_change


_FACE_NEIGHBORS = ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0))
_ZERO = np.float32(1.0e-10)
_LESION_OBJECT_TYPES = {"lesion", "blackhole"}


@dataclass(frozen=True)
class LesionToadsGrowingInitialization:
    """Availability and min-heap entries after ``initGrowingLabels``."""

    available: NDArray[np.bool_]
    queue: ToadsPriorityHeap4D


@dataclass(frozen=True)
class LesionToadsGrowingEvent:
    """Image state after processing one selected source queue item."""

    segmentation: NDArray[np.int16]
    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    ordering: NDArray[np.float32]
    class_volume: NDArray[np.float32]
    popped_item: ToadsHeapItem
    topology_accepted: bool | None
    inserted_neighbors: tuple[ToadsHeapItem, ...]


@dataclass(frozen=True)
class LesionToadsGrowingPass:
    """Final state and selected-source event trace after queue exhaustion."""

    segmentation: NDArray[np.int16]
    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    ordering: NDArray[np.float32]
    class_volume: NDArray[np.float32]
    events: tuple[LesionToadsGrowingEvent, ...]
    remaining_queue: tuple[ToadsHeapItem, ...]


@dataclass(frozen=True)
class LesionToadsThinningInitialization:
    """State and max-heap entries created by ``initThinningLabels``."""

    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    class_volume: NDArray[np.float32]
    boundary_counts: NDArray[np.int32]
    inserted_items: tuple[ToadsHeapItem, ...]
    queue: ToadsPriorityHeap4D


@dataclass(frozen=True)
class LesionToadsThinningTransition:
    """Discrete and volume state after one selected LesionTOADS pop."""

    segmentation_label: int
    skeleton_label: int
    class_volume: NDArray[np.float32]
    available_at_center: NDArray[np.bool_]
    ordering_at_center: np.float32
    candidate_label: int | None
    topology_accepted: bool | None
    processed: bool


@dataclass(frozen=True)
class LesionToadsThinningSchedule:
    """Availability and heap insertions from one thinning post-removal loop."""

    available: NDArray[np.bool_]
    inserted_items: tuple[ToadsHeapItem, ...]


def schedule_lesion_toads_thinning_neighbors(
    memberships: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    available: ArrayLike,
    queue: ToadsPriorityHeap4D,
    removed_voxel: tuple[int, int, int],
    removed_class: int,
    popped_priority: float,
    class_volume: float,
    minimum_class_volume: float,
    limit: float,
    *,
    use_lesion_max: bool = False,
    zero: float = 1.0e-10,
) -> LesionToadsThinningSchedule:
    """Schedule eligible face neighbors using LesionTOADS membership scores.

    This is the selected post-removal loop in ``propagateThinning``. Unlike
    the base TOADS scheduler, LesionTOADS does not special-case a neighbor
    with one same-label adjacent voxel; every available neighbor gets the
    ``val - 1/maxMembershipRatio`` priority using lesion/black-hole scores.
    """
    mems = np.asarray(memberships, dtype=np.float32)
    lesion_mems = np.asarray(lesion_membership, dtype=np.float32)
    black_hole_mems = np.asarray(black_hole_membership, dtype=np.float32)
    state = np.asarray(available)
    types = tuple(object_types)
    if mems.ndim != 4 or mems.shape[-1] < 2 or mems.shape[-1] > 128:
        raise ValueError("memberships must be 4-D with 2 to 128 classes")
    if len(types) != mems.shape[-1]:
        raise ValueError("one object type is required per membership class")
    spatial = mems.shape[:3]
    if any(size < 3 for size in spatial):
        raise ValueError("membership grids must have every dimension at least 3")
    if lesion_mems.shape != spatial or black_hole_mems.shape != spatial:
        raise ValueError("lesion and black-hole membership maps must match memberships")
    if state.shape != (*spatial, mems.shape[-1]) or state.dtype != np.bool_:
        raise ValueError("available must be a boolean grid matching memberships")
    if not all(np.isfinite(array).all() for array in (mems, lesion_mems, black_hole_mems)):
        raise ValueError("membership arrays must be finite")
    if not isinstance(queue, ToadsPriorityHeap4D) or queue.mode != "max":
        raise ValueError("queue must be a TOADS max heap")
    if (not isinstance(removed_class, (int, np.integer))
            or isinstance(removed_class, (bool, np.bool_))
            or not 0 <= int(removed_class) < mems.shape[-1]):
        raise ValueError("removed_class must be a valid membership class")
    if len(removed_voxel) != 3:
        raise ValueError("removed_voxel must contain three coordinates")
    x, y, z = (int(value) for value in removed_voxel)
    if not (1 <= x < spatial[0] - 1 and 1 <= y < spatial[1] - 1 and 1 <= z < spatial[2] - 1):
        raise ValueError("removed_voxel must be interior")
    zero32 = np.float32(zero)
    value32 = np.float32(popped_priority)
    volume32 = np.float32(class_volume)
    minimum32 = np.float32(minimum_class_volume)
    limit32 = np.float32(limit)
    if (not all(np.isfinite(item) for item in (zero32, value32, volume32, minimum32, limit32))
            or zero32 <= np.float32(0.0) or limit32 < np.float32(0.0)):
        raise ValueError("priority/volume/limit must be finite; ZERO must be positive and limit nonnegative")

    updated = state.copy()
    inserted: list[ToadsHeapItem] = []
    # This nested order is the source caller's six face sequence and controls
    # equal-priority insertions into Java's non-stable binary heap.
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                if i * i + j * j + k * k != 1:
                    continue
                xi, yj, zk = x + i, y + j, z + k
                if not updated[xi, yj, zk, removed_class]:
                    continue
                best = _best_membership_at(
                    xi, yj, zk, mems, types, lesion_mems, black_hole_mems,
                    use_lesion_max=bool(use_lesion_max),
                )
                ratio = _membership_ratio_at(
                    xi, yj, zk, int(removed_class), best, mems, types,
                    lesion_mems, black_hole_mems,
                    use_lesion_max=bool(use_lesion_max), zero=zero32,
                )
                order = np.float32(value32 - np.float32(np.float32(1.0) / ratio))
                if volume32 > minimum32 and order > np.float32(-limit32):
                    updated[xi, yj, zk, removed_class] = False
                    queue.add(order, xi, yj, zk, int(removed_class))
                    inserted.append(ToadsHeapItem(order, xi, yj, zk, int(removed_class)))
    return LesionToadsThinningSchedule(updated, tuple(inserted))


def lesion_toads_thinning_voxel_transition(
    segmentation_neighborhood: ArrayLike,
    skeleton_neighborhood: ArrayLike,
    memberships: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    lesion_membership_neighborhood: ArrayLike,
    black_hole_membership_neighborhood: ArrayLike,
    transformed_shape_neighborhood: ArrayLike,
    class_volume: ArrayLike,
    minimum_class_volume: ArrayLike,
    available_at_center: ArrayLike,
    ordering_at_center: float,
    queue_value: float,
    *,
    lesion_class: int,
    processing_class: int,
    use_lesion_max: bool = False,
    empty_label: int = -1,
) -> LesionToadsThinningTransition:
    """Apply one selected LesionTOADS ``propagateThinning`` removal event.

    The caller has already removed the heap head. This function ports the
    selected local target choice, topology gate, volume decrement, and center
    state updates. The shape-prior fallback intentionally follows the
    LesionTOADS source (largest transformed shape value, without TOADS'
    relation-matrix check); queue-neighbor insertion is a separate operation.
    """
    labels = np.asarray(segmentation_neighborhood)
    skeleton = np.asarray(skeleton_neighborhood)
    mems = np.asarray(memberships, dtype=np.float32)
    lesion_mems = np.asarray(lesion_membership_neighborhood, dtype=np.float32)
    black_hole_mems = np.asarray(black_hole_membership_neighborhood, dtype=np.float32)
    shapes = np.asarray(transformed_shape_neighborhood, dtype=np.float32)
    volumes = np.asarray(class_volume, dtype=np.float32)
    minimums = np.asarray(minimum_class_volume, dtype=np.float32)
    available = np.asarray(available_at_center)
    types = tuple(object_types)
    if labels.shape != (3, 3, 3) or skeleton.shape != (3, 3, 3):
        raise ValueError("segmentation and skeleton neighborhoods must have shape (3, 3, 3)")
    if not np.issubdtype(labels.dtype, np.integer) or not np.issubdtype(skeleton.dtype, np.integer):
        raise ValueError("segmentation and skeleton labels must be integers")
    classes = len(types)
    if classes < 2 or classes > 128:
        raise ValueError("object_types must contain between 2 and 128 classes")
    if mems.shape != (classes,) or not np.isfinite(mems).all():
        raise ValueError("memberships must be a finite class vector")
    if lesion_mems.shape != (3, 3, 3) or black_hole_mems.shape != (3, 3, 3):
        raise ValueError("lesion and black-hole neighborhoods must have shape (3, 3, 3)")
    if shapes.shape != (3, 3, 3, classes) or not np.isfinite(shapes).all():
        raise ValueError("transformed shapes must be finite with shape (3, 3, 3, classes)")
    if not np.isfinite(lesion_mems).all() or not np.isfinite(black_hole_mems).all():
        raise ValueError("lesion and black-hole memberships must be finite")
    if volumes.shape != (classes,) or minimums.shape != (classes,):
        raise ValueError("class volume arrays must match object_types")
    if available.shape != (classes,) or available.dtype != np.bool_:
        raise ValueError("available_at_center must be a boolean vector matching object_types")
    if not np.isfinite(volumes).all() or not np.isfinite(minimums).all():
        raise ValueError("class volumes and minimums must be finite")
    for name, index in (("lesion_class", lesion_class), ("processing_class", processing_class)):
        if isinstance(index, (bool, np.bool_)) or not isinstance(index, (int, np.integer)):
            raise ValueError(f"{name} must be an integer class index")
        if not 0 <= int(index) < classes:
            raise ValueError(f"{name} must be a valid class index")
    if not isinstance(use_lesion_max, (bool, np.bool_)):
        raise ValueError("use_lesion_max must be boolean")
    if not np.isfinite([ordering_at_center, queue_value]).all():
        raise ValueError("ordering values must be finite")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation labels must be valid class indices")
    if np.any((skeleton != empty_label) & ((skeleton < 0) | (skeleton >= classes))):
        raise ValueError("skeleton labels must be valid class indices or empty_label")

    center_label = int(labels[1, 1, 1])
    center_skeleton = int(skeleton[1, 1, 1])
    updated_volumes = volumes.copy()
    updated_available = available.copy()
    if center_skeleton != processing_class or not (
        volumes[processing_class] > minimums[processing_class]
    ):
        return LesionToadsThinningTransition(
            center_label, center_skeleton, updated_volumes, updated_available,
            np.float32(ordering_at_center), None, None, False,
        )

    target = _best_membership_at(
        0, 0, 0,
        mems.reshape((1, 1, 1, classes)),
        types,
        lesion_mems[1:2, 1:2, 1:2],
        black_hole_mems[1:2, 1:2, 1:2],
        use_lesion_max=bool(use_lesion_max),
    )
    direct = False
    next_label: int | None = None
    best_prior = np.float32(0.0)
    if target != processing_class:
        for dx, dy, dz in _FACE_NEIGHBORS:
            neighbor = int(labels[1 + dx, 1 + dy, 1 + dz])
            if neighbor == target:
                direct = True
            elif neighbor != processing_class:
                prior = shapes[1 + dx, 1 + dy, 1 + dz, target]
                if prior > best_prior:
                    next_label = neighbor
                    best_prior = prior
    if target == processing_class:
        candidate = None
        topology_accepted = None
    elif direct:
        candidate = target
        topology_accepted = toads_homeomorphic_label_change(labels, target)
    elif next_label is not None:
        candidate = next_label
        topology_accepted = toads_homeomorphic_label_change(labels, next_label)
    else:
        candidate = None
        topology_accepted = None

    updated_available[processing_class] = False
    updated_ordering = np.float32(queue_value)
    if processing_class == lesion_class:
        if use_lesion_max:
            decrement = _lesion_score(
                mems[processing_class], lesion_mems[1, 1, 1], black_hole_mems[1, 1, 1],
                use_lesion_max=True,
            )
        else:
            decrement = np.float32(
                np.float32(mems[processing_class] + lesion_mems[1, 1, 1])
                + black_hole_mems[1, 1, 1]
            )
    else:
        decrement = mems[processing_class]
    updated_volumes[processing_class] = np.float32(
        volumes[processing_class] - decrement
    )
    new_label = (
        int(candidate)
        if candidate is not None and topology_accepted
        else center_label
    )
    return LesionToadsThinningTransition(
        new_label, int(empty_label), updated_volumes, updated_available,
        updated_ordering, candidate, topology_accepted, True,
    )


def _lesion_score(
    ordinary_membership: np.float32,
    lesion_membership: np.float32,
    black_hole_membership: np.float32,
    *,
    use_lesion_max: bool,
) -> np.float32:
    if use_lesion_max:
        # VandyCRUISE Numerics.max(a,b,c) is max(a, max(b,c)), with strict
        # greater-than comparisons; values here are finite and this preserves
        # its float32 result and tie behavior.
        inner = lesion_membership if lesion_membership > black_hole_membership else black_hole_membership
        return ordinary_membership if ordinary_membership > inner else inner
    return np.float32(
        np.float32(ordinary_membership + lesion_membership) + black_hole_membership
    )


def _best_membership_at(
    x: int,
    y: int,
    z: int,
    memberships: NDArray[np.float32],
    object_types: tuple[str, ...],
    lesion_membership: NDArray[np.float32],
    black_hole_membership: NDArray[np.float32],
    *,
    use_lesion_max: bool,
) -> int:
    best = 0
    running_value = np.float32(0.0)
    for candidate, object_type in enumerate(object_types):
        if object_type in _LESION_OBJECT_TYPES:
            score = _lesion_score(
                memberships[x, y, z, candidate],
                lesion_membership[x, y, z],
                black_hole_membership[x, y, z],
                use_lesion_max=use_lesion_max,
            )
        else:
            score = memberships[x, y, z, candidate]
        if score > running_value:
            best = candidate
            # Source quirk: compare the combined lesion score but store only
            # the ordinary membership as the next running maximum.
            running_value = memberships[x, y, z, candidate]
    return best


def _membership_ratio_at(
    x: int,
    y: int,
    z: int,
    label: int,
    best: int,
    memberships: NDArray[np.float32],
    object_types: tuple[str, ...],
    lesion_membership: NDArray[np.float32],
    black_hole_membership: NDArray[np.float32],
    *,
    use_lesion_max: bool,
    zero: np.float32,
) -> np.float32:
    if label == best:
        return np.float32(1.0)
    if object_types[label] in _LESION_OBJECT_TYPES:
        current = _lesion_score(
            memberships[x, y, z, label],
            lesion_membership[x, y, z],
            black_hole_membership[x, y, z],
            use_lesion_max=use_lesion_max,
        )
    else:
        current = memberships[x, y, z, label]
    if object_types[best] in _LESION_OBJECT_TYPES:
        maximum = _lesion_score(
            memberships[x, y, z, best],
            lesion_membership[x, y, z],
            black_hole_membership[x, y, z],
            use_lesion_max=use_lesion_max,
        )
    else:
        maximum = memberships[x, y, z, best]
    numerator = np.float32(maximum + zero)
    denominator = np.float32(current + zero)
    return np.float32(numerator / denominator)


def initialize_lesion_toads_thinning_frontier(
    segmentation: ArrayLike,
    memberships: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    initial_skeleton: ArrayLike,
    initial_available: ArrayLike,
    *,
    lesion_class: int,
    use_lesion_max: bool = False,
    zero: float = 1.0e-10,
) -> LesionToadsThinningInitialization:
    """Initialize the selected LesionTOADS thinning state and max heap.

    This reproduces only ``initThinningLabels`` in
    ``LesionToadSegmentation``. It resets class volumes, updates the interior
    skeleton and per-class availability, then queues six-face boundaries
    with priority ``-1 / maxMembershipRatio``. The selected Java
    ``bestMembership`` quirk and float32 evaluation order are shared with the
    growing-frontier adapter above. Outer planes remain as supplied because
    the source initialization loops only over interior voxels.
    """
    labels = np.asarray(segmentation)
    mems = np.asarray(memberships, dtype=np.float32)
    lesion_mems = np.asarray(lesion_membership, dtype=np.float32)
    black_hole_mems = np.asarray(black_hole_membership, dtype=np.float32)
    skeleton = np.asarray(initial_skeleton)
    available = np.asarray(initial_available)
    types = tuple(object_types)
    if labels.ndim != 3 or any(size < 3 for size in labels.shape):
        raise ValueError("segmentation must be a 3-D grid with every dimension at least 3")
    if not np.issubdtype(labels.dtype, np.integer) or np.issubdtype(labels.dtype, np.bool_):
        raise ValueError("segmentation labels must be integers")
    if mems.ndim != 4 or mems.shape[:3] != labels.shape or mems.shape[-1] < 1:
        raise ValueError("memberships must match segmentation and contain a class axis")
    classes = mems.shape[-1]
    if classes > 128 or len(types) != classes:
        raise ValueError("one object type is required per Java signed-byte class")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation contains a class index outside the membership axis")
    if lesion_mems.shape != labels.shape or black_hole_mems.shape != labels.shape:
        raise ValueError("lesion and black-hole maps must match segmentation")
    if skeleton.shape != labels.shape or not np.issubdtype(skeleton.dtype, np.integer):
        raise ValueError("initial_skeleton must be an integer grid matching segmentation")
    if available.shape != (*labels.shape, classes) or available.dtype != np.bool_:
        raise ValueError("initial_available must be a boolean grid matching segmentation and classes")
    if not all(np.isfinite(array).all() for array in (mems, lesion_mems, black_hole_mems)):
        raise ValueError("membership arrays must contain only finite values")
    if isinstance(lesion_class, (bool, np.bool_)) or not isinstance(lesion_class, (int, np.integer)):
        raise ValueError("lesion_class must be an integer class index")
    if not 0 <= int(lesion_class) < classes:
        raise ValueError("lesion_class is outside the membership class axis")
    if not isinstance(use_lesion_max, (bool, np.bool_)):
        raise ValueError("use_lesion_max must be boolean")
    zero_value = np.float32(zero)
    if not np.isfinite(zero_value) or zero_value <= np.float32(0.0):
        raise ValueError("zero must be finite and positive")

    skeleton_out = skeleton.astype(np.int16, copy=True)
    available_out = available.copy()
    volume = np.zeros(classes, dtype=np.float32)
    counts = np.zeros(classes, dtype=np.int32)
    nx, ny, nz = labels.shape
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                current = int(labels[x, y, z])
                skeleton_out[x, y, z] = current
                available_out[x, y, z, :] = False
                available_out[x, y, z, current] = True
                if types[current] in _LESION_OBJECT_TYPES:
                    if use_lesion_max:
                        increment = _lesion_score(
                            mems[x, y, z, current],
                            lesion_mems[x, y, z],
                            black_hole_mems[x, y, z],
                            use_lesion_max=True,
                        )
                    else:
                        increment = np.float32(
                            np.float32(mems[x, y, z, current] + lesion_mems[x, y, z])
                            + black_hole_mems[x, y, z]
                        )
                else:
                    increment = mems[x, y, z, current]
                volume[current] = np.float32(volume[current] + increment)

    heap = ToadsPriorityHeap4D(max(1, int(np.prod(labels.shape)) * classes), mode="max")
    inserted: list[ToadsHeapItem] = []
    face_neighbors = ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1))
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                current = int(labels[x, y, z])
                if not any(
                    labels[x + dx, y + dy, z + dz] != current
                    for dx, dy, dz in face_neighbors
                ):
                    continue
                available_out[x, y, z, current] = False
                best = _best_membership_at(
                    x,
                    y,
                    z,
                    mems,
                    types,
                    lesion_mems,
                    black_hole_mems,
                    use_lesion_max=bool(use_lesion_max),
                )
                ratio = _membership_ratio_at(
                    x,
                    y,
                    z,
                    current,
                    best,
                    mems,
                    types,
                    lesion_mems,
                    black_hole_mems,
                    use_lesion_max=bool(use_lesion_max),
                    zero=zero_value,
                )
                priority = np.float32(-np.float32(np.float32(1.0) / ratio))
                item = ToadsHeapItem(priority, x, y, z, current)
                heap.add(priority, x, y, z, current)
                inserted.append(item)
                counts[current] += 1

    return LesionToadsThinningInitialization(
        skeleton=skeleton_out,
        available=available_out,
        class_volume=volume,
        boundary_counts=counts,
        inserted_items=tuple(inserted),
        queue=heap,
    )


def initialize_lesion_toads_growing_frontier(
    skeleton: ArrayLike,
    memberships: ArrayLike,
    ordering: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    *,
    initial_available: ArrayLike | None = None,
    use_lesion_max: bool = False,
    empty_label: int = -1,
    zero: float = 1.0e-10,
) -> LesionToadsGrowingInitialization:
    """Seed the selected LesionTOADS growing queue from a class skeleton.

    Inputs use ``(x, y, z)`` voxels and ``(x, y, z, class)`` memberships.
    Interior voxels are reset available for every class iff their skeleton
    label is empty. Each available voxel that touches class ``k`` on a
    six-face is removed from availability for ``k`` and inserted with priority
    ``ordering + maxMembershipRatio``. Border availability is preserved from
    the input but must be false for the source's unguarded neighbor reads.

    The selected Java ``bestMembership`` has a source quirk: it compares the
    adjusted lesion score against its running maximum, but when it wins it
    stores the ordinary membership (not that adjusted score) as the next
    maximum. This function retains that operation order. Only source object
    types ``lesion`` and ``blackhole`` receive the special score; optimized
    lesion/black-hole types do not in these two methods.
    """
    skel = np.asarray(skeleton)
    mems = np.asarray(memberships, dtype=np.float32)
    orderings = np.asarray(ordering, dtype=np.float32)
    lesion_mems = np.asarray(lesion_membership, dtype=np.float32)
    black_hole_mems = np.asarray(black_hole_membership, dtype=np.float32)
    types = tuple(object_types)

    if skel.ndim != 3 or any(size < 3 for size in skel.shape):
        raise ValueError("skeleton must be a 3-D grid with every dimension at least 3")
    if not np.issubdtype(skel.dtype, np.integer):
        raise ValueError("skeleton labels must be integers")
    if mems.ndim != 4 or mems.shape[:3] != skel.shape or mems.shape[-1] < 2:
        raise ValueError("memberships must match skeleton and contain at least two classes")
    classes = int(mems.shape[-1])
    if classes > 127:
        raise ValueError("class indices must fit the Java signed-byte loop")
    if len(types) != classes or not all(isinstance(value, str) and value for value in types):
        raise ValueError("one nonempty object type is required for every class")
    if orderings.shape != skel.shape:
        raise ValueError("ordering must match the skeleton grid")
    for name, values in (
        ("memberships", mems),
        ("ordering", orderings),
        ("lesion_membership", lesion_mems),
        ("black_hole_membership", black_hole_mems),
    ):
        if not np.isfinite(values).all():
            raise ValueError(f"{name} must contain only finite values")
    if lesion_mems.shape != skel.shape or black_hole_mems.shape != skel.shape:
        raise ValueError("lesion and black-hole memberships must match skeleton")
    if not isinstance(empty_label, (int, np.integer)):
        raise ValueError("empty_label must be an integer")
    if np.any((skel != empty_label) & ((skel < 0) | (skel >= classes))):
        raise ValueError("skeleton labels must be valid class indices or empty_label")
    zero_value = np.float32(zero)
    if not np.isfinite(zero_value) or zero_value < np.float32(0.0):
        raise ValueError("zero must be finite and nonnegative")

    if initial_available is None:
        available = np.zeros(skel.shape + (classes,), dtype=np.bool_)
    else:
        supplied = np.asarray(initial_available)
        if supplied.shape != skel.shape + (classes,) or supplied.dtype != np.bool_:
            raise ValueError("initial_available must be boolean and match skeleton and classes")
        available = supplied.copy()
    border = np.zeros(skel.shape, dtype=np.bool_)
    for axis in range(3):
        edge = [slice(None)] * 3
        edge[axis] = 0
        border[tuple(edge)] = True
        edge[axis] = -1
        border[tuple(edge)] = True
    if np.any(available[border]):
        raise ValueError("initial border availability must be false for source-safe neighbor reads")

    nx, ny, nz = skel.shape
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                available[x, y, z, :] = int(skel[x, y, z]) == empty_label

    queue = ToadsPriorityHeap4D(nx * ny + ny * nz + nz * nx, "min")
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for label in range(classes):
                    if not available[x, y, z, label]:
                        continue
                    if not any(
                        int(skel[x + dx, y + dy, z + dz]) == label
                        for dx, dy, dz in _FACE_NEIGHBORS
                    ):
                        continue
                    available[x, y, z, label] = False

                    best = _best_membership_at(
                        x,
                        y,
                        z,
                        mems,
                        types,
                        lesion_mems,
                        black_hole_mems,
                        use_lesion_max=use_lesion_max,
                    )
                    ratio = _membership_ratio_at(
                        x,
                        y,
                        z,
                        label,
                        best,
                        mems,
                        types,
                        lesion_mems,
                        black_hole_mems,
                        use_lesion_max=use_lesion_max,
                        zero=zero_value,
                    )
                    priority = np.float32(orderings[x, y, z] + ratio)
                    queue.add(priority, x, y, z, label)

    return LesionToadsGrowingInitialization(available=available, queue=queue)


def process_lesion_toads_growing_event(
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    memberships: ArrayLike,
    available: ArrayLike,
    ordering: ArrayLike,
    class_volume: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    queue: ToadsPriorityHeap4D,
    *,
    lesion_class: int,
    use_lesion_max: bool = False,
    empty_label: int = -1,
    zero: float = 1.0e-10,
) -> LesionToadsGrowingEvent:
    """Apply one selected ``LesionToadSegmentation.propagateGrowing`` pop.

    The event composes the already checked 18/6 local topology predicate and
    P39q lesion-aware scoring with the source center updates and six-face
    neighbor scheduling. It stops after one pop. A rejected topology event
    re-enables that center/class availability; a stale event whose skeleton
    is already occupied is ignored without invoking topology.
    """
    labels = np.asarray(segmentation)
    skel = np.asarray(skeleton)
    mems = np.asarray(memberships, dtype=np.float32)
    avail = np.asarray(available)
    orderings = np.asarray(ordering, dtype=np.float32)
    volumes = np.asarray(class_volume, dtype=np.float32)
    lesion_mems = np.asarray(lesion_membership, dtype=np.float32)
    black_hole_mems = np.asarray(black_hole_membership, dtype=np.float32)
    types = tuple(object_types)
    if labels.ndim != 3 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be a 3-D integer grid")
    if any(size < 3 for size in labels.shape):
        raise ValueError("segmentation dimensions must each be at least 3")
    if skel.shape != labels.shape or not np.issubdtype(skel.dtype, np.integer):
        raise ValueError("skeleton must be an integer grid matching segmentation")
    if mems.ndim != 4 or mems.shape[:3] != labels.shape or mems.shape[-1] < 2:
        raise ValueError("memberships must be 4-D and match segmentation")
    classes = int(mems.shape[-1])
    if classes > 127:
        raise ValueError("class indices must fit the Java signed-byte loop")
    if len(types) != classes or not all(isinstance(value, str) and value for value in types):
        raise ValueError("one nonempty object type is required for every class")
    if avail.shape != labels.shape + (classes,) or avail.dtype != np.bool_:
        raise ValueError("available must be boolean and match segmentation and classes")
    if orderings.shape != labels.shape or volumes.shape != (classes,):
        raise ValueError("ordering and class volumes must match the image and classes")
    if lesion_mems.shape != labels.shape or black_hole_mems.shape != labels.shape:
        raise ValueError("lesion and black-hole memberships must match segmentation")
    for name, values in (
        ("memberships", mems),
        ("ordering", orderings),
        ("class_volume", volumes),
        ("lesion_membership", lesion_mems),
        ("black_hole_membership", black_hole_mems),
    ):
        if not np.isfinite(values).all():
            raise ValueError(f"{name} must contain only finite values")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation labels must be valid class indices")
    if np.any((skel != empty_label) & ((skel < 0) | (skel >= classes))):
        raise ValueError("skeleton labels must be valid class indices or empty_label")
    if not isinstance(lesion_class, (int, np.integer)) or not -1 <= lesion_class < classes:
        raise ValueError("lesion_class must be -1 or a valid class index")
    zero_value = np.float32(zero)
    if not np.isfinite(zero_value) or zero_value < np.float32(0.0):
        raise ValueError("zero must be finite and nonnegative")
    if not isinstance(queue, ToadsPriorityHeap4D) or queue.mode != "min":
        raise ValueError("queue must be a TOADS min heap")
    if not queue:
        raise ValueError("cannot process an empty LesionTOADS growing queue")

    item = queue.peek_first()
    x, y, z, label = item.x, item.y, item.z, item.label
    if not (
        1 <= x < labels.shape[0] - 1
        and 1 <= y < labels.shape[1] - 1
        and 1 <= z < labels.shape[2] - 1
        and 0 <= label < classes
    ):
        raise ValueError("queue item must have an interior coordinate and valid class")
    border = np.zeros(labels.shape, dtype=np.bool_)
    for axis in range(3):
        edge = [slice(None)] * 3
        edge[axis] = 0
        border[tuple(edge)] = True
        edge[axis] = -1
        border[tuple(edge)] = True
    if np.any(avail[border]):
        raise ValueError("border availability must be false for safe source neighbor updates")

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
            if label == lesion_class:
                membership_increment = _lesion_score(
                    mems[x, y, z, label],
                    lesion_mems[x, y, z],
                    black_hole_mems[x, y, z],
                    use_lesion_max=use_lesion_max,
                )
            else:
                membership_increment = mems[x, y, z, label]
            updated_volumes[label] = np.float32(
                updated_volumes[label] + membership_increment
            )
            updated_labels[x, y, z] = np.int16(label)

            for dx, dy, dz in _FACE_NEIGHBORS:
                xi, yj, zl = x + dx, y + dy, z + dz
                if not updated_available[xi, yj, zl, label]:
                    continue
                best = _best_membership_at(
                    xi,
                    yj,
                    zl,
                    mems,
                    types,
                    lesion_mems,
                    black_hole_mems,
                    use_lesion_max=use_lesion_max,
                )
                ratio = _membership_ratio_at(
                    xi,
                    yj,
                    zl,
                    label,
                    best,
                    mems,
                    types,
                    lesion_mems,
                    black_hole_mems,
                    use_lesion_max=use_lesion_max,
                    zero=zero_value,
                )
                priority = np.float32(popped.value + ratio)
                updated_available[xi, yj, zl, label] = False
                queue.add(priority, xi, yj, zl, label)
                inserted.append(ToadsHeapItem(priority, xi, yj, zl, label))
        else:
            updated_available[x, y, z, label] = True

    return LesionToadsGrowingEvent(
        segmentation=updated_labels,
        skeleton=updated_skeleton,
        available=updated_available,
        ordering=updated_ordering,
        class_volume=updated_volumes,
        popped_item=popped,
        topology_accepted=topology_accepted,
        inserted_neighbors=tuple(inserted),
    )


def propagate_lesion_toads_growing(
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    memberships: ArrayLike,
    ordering: ArrayLike,
    class_volume: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    *,
    lesion_class: int,
    use_lesion_max: bool = False,
    empty_label: int = -1,
    zero: float = 1.0e-10,
) -> LesionToadsGrowingPass:
    """Run the selected lesion-aware growing pass until its min-heap empties.

    This composes the separately checked ``initGrowingLabels`` and one-pop
    ``propagateGrowing`` operations. Returned arrays do not alias the caller's
    arrays. The event trace is useful for source comparison and provenance.
    No lesion-membership estimation, competition, or thinning is performed.
    """
    labels = np.asarray(segmentation)
    skel = np.asarray(skeleton)
    mems = np.asarray(memberships, dtype=np.float32)
    orderings = np.asarray(ordering, dtype=np.float32)
    volumes = np.asarray(class_volume, dtype=np.float32)
    lesion_mems = np.asarray(lesion_membership, dtype=np.float32)
    black_hole_mems = np.asarray(black_hole_membership, dtype=np.float32)
    types = tuple(object_types)
    initialized = initialize_lesion_toads_growing_frontier(
        skel,
        mems,
        orderings,
        types,
        lesion_mems,
        black_hole_mems,
        use_lesion_max=use_lesion_max,
        empty_label=empty_label,
        zero=zero,
    )
    current_labels = labels.astype(np.int16, copy=True)
    current_skeleton = skel.astype(np.int16, copy=True)
    current_available = initialized.available
    current_ordering = orderings.copy()
    current_volume = volumes.copy()
    events: list[LesionToadsGrowingEvent] = []
    while initialized.queue:
        event = process_lesion_toads_growing_event(
            current_labels,
            current_skeleton,
            mems,
            current_available,
            current_ordering,
            current_volume,
            types,
            lesion_mems,
            black_hole_mems,
            initialized.queue,
            lesion_class=lesion_class,
            use_lesion_max=use_lesion_max,
            empty_label=empty_label,
            zero=zero,
        )
        events.append(event)
        current_labels = event.segmentation
        current_skeleton = event.skeleton
        current_available = event.available
        current_ordering = event.ordering
        current_volume = event.class_volume

    return LesionToadsGrowingPass(
        segmentation=current_labels,
        skeleton=current_skeleton,
        available=current_available,
        ordering=current_ordering,
        class_volume=current_volume,
        events=tuple(events),
        remaining_queue=tuple(),
    )
