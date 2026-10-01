"""Selected LesionTOADS narrow-band competition queue pass.

This is a source-scoped port of ``initCompetitionLabels`` and
``propagateCompetition`` from VandyCRUISE ``LesionToadSegmentation``. It
consumes initialized segmentation, membership, lesion/black-hole membership,
and class-volume state; it does not estimate those inputs or run the wrapper's
statistics, registration, growth, or stopping stages.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.toads_heap import ToadsHeapItem, ToadsPriorityHeap4D
from pymacruise.toads_topology import toads_homeomorphic_label_change


_ZERO = np.float32(1e-10)
_FACE_OFFSETS = (
    (-1, 0, 0), (0, -1, 0), (0, 0, -1),
    (0, 0, 1), (0, 1, 0), (1, 0, 0),
)


@dataclass(frozen=True)
class LesionToadsCompetitionInitialization:
    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    inserted_items: tuple[ToadsHeapItem, ...]
    queue: ToadsPriorityHeap4D


@dataclass(frozen=True)
class LesionToadsCompetitionEvent:
    popped_item: ToadsHeapItem
    processed: bool
    candidate_label: int
    topology_accepted: bool | None
    inserted_neighbors: tuple[ToadsHeapItem, ...]


@dataclass(frozen=True)
class LesionToadsCompetitionPass:
    segmentation: NDArray[np.int16]
    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    ordering: NDArray[np.float32]
    class_volume: NDArray[np.float32]
    initialization: LesionToadsCompetitionInitialization
    events: tuple[LesionToadsCompetitionEvent, ...]
    remaining_queue: tuple[ToadsHeapItem, ...]


def initialize_lesion_toads_competition(
    segmentation: ArrayLike,
    memberships: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    *,
    lesion_class: int,
    object_types: tuple[str, ...],
) -> LesionToadsCompetitionInitialization:
    """Create the selected max-priority transition queue and reset skeleton.

    Arrays use x/y/z[/class] order. ``object_types`` entries are ``lesion``,
    ``blackhole``, or any other selected class type. The returned state owns
    new skeleton/availability arrays; caller arrays are not modified.
    """
    labels, mems, lesion, blackhole, classes = _validate_inputs(
        segmentation, memberships, lesion_membership, black_hole_membership,
        lesion_class, object_types,
    )
    shape = labels.shape
    available = np.ones((*shape, classes), dtype=np.bool_)
    skeleton = labels.copy()
    queue = ToadsPriorityHeap4D(max(1, int(np.prod(shape)) * 6), mode="max")
    inserted: list[ToadsHeapItem] = []

    nx, ny, nz = shape
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                for k in range(classes):
                    if int(labels[x, y, z]) != k:
                        continue
                    for dx, dy, dz in _FACE_OFFSETS:
                        xi, yj, zl = x + dx, y + dy, z + dz
                        if not (0 <= xi < nx and 0 <= yj < ny and 0 <= zl < nz):
                            continue
                        m = int(labels[xi, yj, zl])
                        if m == k:
                            continue
                        order = _competition_initial_order(
                            mems[xi, yj, zl], lesion[xi, yj, zl],
                            blackhole[xi, yj, zl], k, m, lesion_class,
                        )
                        if order > 0:
                            packed = k + classes * m
                            queue.add(order, x, y, z, packed)
                            inserted.append(ToadsHeapItem(
                                np.float32(order), x, y, z, packed
                            ))
    return LesionToadsCompetitionInitialization(
        skeleton, available, tuple(inserted), queue
    )


def propagate_lesion_toads_competition(
    segmentation: ArrayLike,
    memberships: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    ordering: ArrayLike,
    class_volume: ArrayLike,
    *,
    lesion_class: int,
    object_types: tuple[str, ...],
) -> LesionToadsCompetitionPass:
    """Run one selected competition pass to queue exhaustion.

    The volume update intentionally preserves the selected Java expression
    ``volume[k] -= mems[k] - lesionMembership`` for a lesion source class.
    """
    labels, mems, lesion, blackhole, classes = _validate_inputs(
        segmentation, memberships, lesion_membership, black_hole_membership,
        lesion_class, object_types,
    )
    source_order = np.asarray(ordering)
    volumes = np.asarray(class_volume)
    if source_order.shape != labels.shape or source_order.dtype != np.float32:
        raise ValueError("ordering must be float32 with segmentation shape")
    if volumes.shape != (classes,) or volumes.dtype != np.float32:
        raise ValueError("class_volume must be float32 with one value per class")
    before = (labels.copy(), mems.copy(), lesion.copy(), blackhole.copy(), source_order.copy(), volumes.copy())

    initialization = initialize_lesion_toads_competition(
        labels, mems, lesion, blackhole, lesion_class=lesion_class,
        object_types=object_types,
    )
    labels = labels.copy()
    skeleton = initialization.skeleton.copy()
    available = initialization.available.copy()
    source_order = source_order.copy()
    volumes = volumes.copy()
    queue = initialization.queue
    events: list[LesionToadsCompetitionEvent] = []
    nx, ny, nz = labels.shape
    # The selected Java method initializes order once, outside the queue loop;
    # special-type branches can intentionally leave its previous value intact.
    neighbor_order_value = np.float32(0.0)

    while queue:
        item = queue.pop_first()
        packed = int(item.label)
        target = packed // classes
        old = packed - target * classes
        inserted: list[ToadsHeapItem] = []
        accepted: bool | None = None
        processed = int(labels[item.x, item.y, item.z]) == old
        if processed:
            neighborhood = labels[
                item.x - 1:item.x + 2,
                item.y - 1:item.y + 2,
                item.z - 1:item.z + 2,
            ]
            accepted = toads_homeomorphic_label_change(neighborhood, target)
            if accepted:
                labels[item.x, item.y, item.z] = target
                skeleton[item.x, item.y, item.z] = target
                source_order[item.x, item.y, item.z] = item.value
                _update_class_volumes(
                    volumes, mems[item.x, item.y, item.z],
                    lesion[item.x, item.y, item.z],
                    blackhole[item.x, item.y, item.z], old, target,
                    lesion_class, object_types,
                )
                for dx, dy, dz in _FACE_OFFSETS:
                    xi, yj, zl = item.x + dx, item.y + dy, item.z + dz
                    if int(labels[xi, yj, zl]) != old:
                        continue
                    computed_order = _competition_neighbor_order(
                        mems[xi, yj, zl], lesion[item.x, item.y, item.z],
                        blackhole[item.x, item.y, item.z], old, target,
                        lesion_class, object_types,
                    )
                    if computed_order is not None:
                        neighbor_order_value = computed_order
                    if neighbor_order_value > 0:
                        queue.add(neighbor_order_value, xi, yj, zl, old)
                        inserted.append(ToadsHeapItem(
                            np.float32(neighbor_order_value), xi, yj, zl, old
                        ))
        events.append(LesionToadsCompetitionEvent(
            item, processed, target, accepted, tuple(inserted)
        ))

    # The source marks every voxel/class available at initialization and does
    # not change availability during competition.
    assert np.array_equal(before[0], np.asarray(segmentation))
    assert np.array_equal(before[1], np.asarray(memberships))
    assert np.array_equal(before[2], np.asarray(lesion_membership))
    assert np.array_equal(before[3], np.asarray(black_hole_membership))
    assert np.array_equal(before[4], np.asarray(ordering))
    assert np.array_equal(before[5], np.asarray(class_volume))
    return LesionToadsCompetitionPass(
        labels, skeleton, available, source_order, volumes, initialization,
        tuple(events), tuple(),
    )


def _validate_inputs(segmentation, memberships, lesion, blackhole, lesion_class, object_types):
    labels = np.asarray(segmentation)
    mems = np.asarray(memberships)
    lesion_values = np.asarray(lesion)
    blackhole_values = np.asarray(blackhole)
    if labels.ndim != 3 or labels.dtype != np.int16:
        raise ValueError("segmentation must be a 3D int16 x/y/z array")
    if mems.ndim != 4 or mems.shape[:3] != labels.shape or mems.dtype != np.float32:
        raise ValueError("memberships must be float32 x/y/z/class and match segmentation")
    classes = mems.shape[-1]
    if lesion_values.shape != labels.shape or lesion_values.dtype != np.float32:
        raise ValueError("lesion_membership must be float32 with segmentation shape")
    if blackhole_values.shape != labels.shape or blackhole_values.dtype != np.float32:
        raise ValueError("black_hole_membership must be float32 with segmentation shape")
    if len(object_types) != classes:
        raise ValueError("object_types must contain one entry per membership class")
    if not 0 <= lesion_class < classes:
        raise ValueError("lesion_class must index the membership classes")
    if any(int(v) < 0 or int(v) >= classes for v in np.unique(labels)):
        raise ValueError("segmentation labels must index the membership classes")
    return labels, mems, lesion_values, blackhole_values, classes


def _competition_initial_order(mems, lesion, blackhole, old, target, lesion_class):
    if target == lesion_class:
        numerator = np.float32(np.float32(mems[target] + lesion) + blackhole)
        denominator = np.float32(mems[old] + _ZERO)
    elif old == lesion_class:
        numerator = np.float32(mems[target])
        denominator = np.float32(np.float32(np.float32(mems[old] + lesion) + blackhole) + _ZERO)
    else:
        numerator = np.float32(mems[target])
        denominator = np.float32(mems[old] + _ZERO)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.float32(numerator / denominator)


def _competition_neighbor_order(mems, lesion, blackhole, old, target, lesion_class, object_types):
    old_type = object_types[old].lower()
    if old == lesion_class:
        if old_type == "lesion":
            denominator = np.float32(np.float32(mems[old] + lesion) + _ZERO)
        elif old_type == "blackhole":
            denominator = np.float32(np.float32(np.float32(mems[old] + lesion) + blackhole) + _ZERO)
        else:
            return None
        numerator = np.float32(mems[target])
    elif target == lesion_class:
        if old_type == "lesion":
            numerator = np.float32(mems[target] + lesion)
        elif old_type == "blackhole":
            numerator = np.float32(np.float32(mems[target] + lesion) + blackhole)
        else:
            return None
        denominator = np.float32(mems[old] + _ZERO)
    else:
        numerator = np.float32(mems[target])
        denominator = np.float32(mems[old] + _ZERO)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.float32(numerator / denominator)


def _update_class_volumes(volumes, mems, lesion, blackhole, old, target, lesion_class, object_types):
    if target == lesion_class:
        increment = np.float32(mems[target] + lesion)
        volumes[target] = np.float32(volumes[target] + increment)
        if object_types[target].lower() == "blackhole":
            volumes[target] = np.float32(volumes[target] + blackhole)
    else:
        volumes[target] = np.float32(volumes[target] + mems[target])
    if old == lesion_class:
        # Intentional source transcription: Java subtracts mems - lesion.
        decrement = np.float32(mems[old] - lesion)
        volumes[old] = np.float32(volumes[old] - decrement)
        if object_types[old].lower() == "blackhole":
            volumes[old] = np.float32(volumes[old] - blackhole)
    else:
        volumes[old] = np.float32(volumes[old] - mems[old])
