"""Sparse-membership growing queue for the selected TOADS approximate branch.

This ports the ``useApprox`` branches in
``MultipleToadSegmentation.initGrowingLabels`` and
``propagateGrowing``. The exact growing implementation remains the reference
for topology checks and source queue semantics; this adapter keeps only the
configured sparse membership slots and never materializes a dense
class-membership volume.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from numpy.typing import ArrayLike

from .toads_heap import ToadsHeapItem, ToadsPriorityHeap4D
from .toads_thinning import (
    ToadsGrowingPass,
    ToadsGrowingQueueInitialization,
    ToadsGrowingQueueStep,
    _FACE_NEIGHBORS,
)
from .toads_topology import toads_homeomorphic_label_change


def _approximate_ratio(
    memberships: np.ndarray,
    best: np.ndarray,
    coordinate: tuple[int, int, int],
    label: int,
    zero: np.float32,
) -> np.float32:
    """Return Java ``maxApproxMembershipRatio`` with float32 operation order."""
    x, y, z = coordinate
    if label == int(best[x, y, z, 0]):
        return np.float32(1.0)
    current = np.float32(0.0)
    for slot in range(1, memberships.shape[-1]):
        if int(best[x, y, z, slot]) == label:
            current = memberships[x, y, z, slot]
    return np.float32(
        np.float32(memberships[x, y, z, 0] + zero)
        / np.float32(current + zero)
    )


def initialize_toads_approximate_growing_queue(
    skeleton: ArrayLike,
    sparse_memberships: ArrayLike,
    best_class_indices: ArrayLike,
    ordering: ArrayLike,
    *,
    class_count: int | None = None,
    initial_available: ArrayLike | None = None,
    empty_label: int = -1,
    zero_threshold: float = 1e-10,
) -> ToadsGrowingQueueInitialization:
    """Initialize selected sparse ``useApprox`` growing candidates."""
    skel = np.asarray(skeleton)
    mems = np.asarray(sparse_memberships, dtype=np.float32)
    best = np.asarray(best_class_indices)
    orderings = np.asarray(ordering, dtype=np.float32)
    if skel.ndim != 3 or any(size < 3 for size in skel.shape) or not np.issubdtype(skel.dtype, np.integer):
        raise ValueError("skeleton must be an integer 3-D grid with every axis at least three")
    if mems.ndim != 4 or mems.shape[:3] != skel.shape or mems.shape[-1] < 2:
        raise ValueError("sparse_memberships must match skeleton and contain at least two slots")
    if best.shape != mems.shape or not np.issubdtype(best.dtype, np.integer):
        raise ValueError("best_class_indices must be an integer grid matching sparse_memberships")
    classes = int(class_count) if class_count is not None else int(np.max(best)) + 1
    if classes < 2 or np.any(best < 0) or np.any(best >= classes):
        raise ValueError("sparse class indices must be valid for at least two classes")
    if not np.isfinite(mems).all():
        raise ValueError("sparse_memberships must be finite")
    if orderings.shape != skel.shape or not np.isfinite(orderings).all():
        raise ValueError("ordering must be finite and match skeleton")
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
    border = np.zeros(skel.shape, dtype=bool)
    for axis in range(3):
        low = [slice(None)] * 3
        high = [slice(None)] * 3
        low[axis], high[axis] = 0, -1
        border[tuple(low)] = True
        border[tuple(high)] = True
    if np.any(available[border]):
        raise ValueError("initial border availability must be false for source-safe face-neighbor reads")

    nx, ny, nz = skel.shape
    updated_skeleton = skel.astype(np.int16, copy=True)
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
                        int(updated_skeleton[x + dx, y + dy, z + dz]) == label
                        for dx, dy, dz in _FACE_NEIGHBORS
                    ):
                        continue
                    available[x, y, z, label] = False
                    ratio = _approximate_ratio(mems, best, (x, y, z), label, zero)
                    priority = np.float32(orderings[x, y, z] + ratio)
                    queue.add(priority, x, y, z, label)
    return ToadsGrowingQueueInitialization(updated_skeleton, available, queue)


def toads_approximate_propagate_growing(
    segmentation: ArrayLike,
    skeleton: ArrayLike,
    sparse_memberships: ArrayLike,
    best_class_indices: ArrayLike,
    ordering: ArrayLike,
    class_volumes: ArrayLike,
    *,
    class_count: int | None = None,
    initial_available: ArrayLike | None = None,
    empty_label: int = -1,
    zero_threshold: float = 1e-10,
    record_events: bool = True,
    progress_callback: Callable[[str], None] | None = None,
    progress_interval: int = 100_000,
) -> ToadsGrowingPass:
    """Initialize and exhaust one sparse approximate TOADS growing pass.

    Sparse volume accumulation follows Java's slot scan and float32 updates.
    On accepted growth, the new voxel's memberships for the grown class are
    added; face neighbors use ``maxApproxMembershipRatio`` and the source
    min-heap. ``record_events=False`` bounds memory for larger images.
    """
    labels_in = np.asarray(segmentation)
    source_skeleton = np.asarray(skeleton)
    mems = np.asarray(sparse_memberships, dtype=np.float32)
    best = np.asarray(best_class_indices)
    orderings = np.asarray(ordering, dtype=np.float32)
    volumes_in = np.asarray(class_volumes, dtype=np.float32)
    if labels_in.ndim != 3 or not np.issubdtype(labels_in.dtype, np.integer):
        raise ValueError("segmentation must be a 3-D integer grid")
    if source_skeleton.shape != labels_in.shape:
        raise ValueError("skeleton must match segmentation")
    if mems.ndim != 4 or mems.shape[:3] != labels_in.shape:
        raise ValueError("sparse_memberships must be 4-D and match segmentation")
    if best.shape != mems.shape or not np.issubdtype(best.dtype, np.integer):
        raise ValueError("best_class_indices must match sparse_memberships")
    classes = int(class_count) if class_count is not None else int(np.max(best)) + 1
    if volumes_in.shape != (classes,) or not np.isfinite(volumes_in).all():
        raise ValueError("class_volumes must be finite and match class_count")
    if orderings.shape != labels_in.shape or not np.isfinite(orderings).all():
        raise ValueError("ordering must be finite and match segmentation")
    if isinstance(progress_interval, bool) or not isinstance(progress_interval, int) or progress_interval < 1:
        raise ValueError("progress_interval must be a positive integer")
    zero = np.float32(zero_threshold)
    if not np.isfinite(zero) or zero < np.float32(0.0):
        raise ValueError("zero_threshold must be finite and nonnegative")
    if not isinstance(empty_label, (int, np.integer)):
        raise ValueError("empty_label must be an integer")

    initialized = initialize_toads_approximate_growing_queue(
        source_skeleton,
        mems,
        best,
        orderings,
        class_count=classes,
        initial_available=initial_available,
        empty_label=empty_label,
        zero_threshold=float(zero),
    )
    current_labels = labels_in.astype(np.int16, copy=True)
    current_skeleton = initialized.skeleton
    available = initialized.available
    current_ordering = orderings.copy()
    volumes = volumes_in.copy()
    queue = initialized.queue
    events: list[ToadsGrowingQueueStep] = []
    processed = 0
    if progress_callback is not None:
        progress_callback(f"approximate growing queue started with {len(queue)} events")

    while queue:
        popped = queue.pop_first()
        x, y, z, label = popped.x, popped.y, popped.z, popped.label
        inserted: list[ToadsHeapItem] = []
        topology_accepted: bool | None = None
        if int(current_skeleton[x, y, z]) == empty_label:
            neighborhood = current_labels[x - 1:x + 2, y - 1:y + 2, z - 1:z + 2]
            topology_accepted = toads_homeomorphic_label_change(neighborhood, label)
            if topology_accepted:
                current_skeleton[x, y, z] = np.int16(label)
                available[x, y, z, :] = False
                current_ordering[x, y, z] = popped.value
                for slot in range(mems.shape[-1]):
                    if int(best[x, y, z, slot]) == label:
                        volumes[label] = np.float32(volumes[label] + mems[x, y, z, slot])
                current_labels[x, y, z] = np.int16(label)
                for dx, dy, dz in _FACE_NEIGHBORS:
                    xi, yj, zl = x + dx, y + dy, z + dz
                    if not available[xi, yj, zl, label]:
                        continue
                    ratio = _approximate_ratio(mems, best, (xi, yj, zl), label, zero)
                    priority = np.float32(popped.value + ratio)
                    available[xi, yj, zl, label] = False
                    queue.add(priority, xi, yj, zl, label)
                    inserted.append(ToadsHeapItem(priority, xi, yj, zl, label))
            else:
                available[x, y, z, label] = True

        if record_events:
            events.append(ToadsGrowingQueueStep(
                current_labels.copy(), current_skeleton.copy(), available.copy(),
                current_ordering.copy(), volumes.copy(), popped, topology_accepted,
                tuple(inserted),
            ))
        processed += 1
        if progress_callback is not None and processed % progress_interval == 0:
            progress_callback(f"approximate growing queue processed {processed} events; {len(queue)} remain")

    if progress_callback is not None:
        progress_callback(f"approximate growing queue complete after {processed} events")
    return ToadsGrowingPass(
        current_labels, current_skeleton, available, current_ordering, volumes,
        tuple(events), (),
    )
