"""Reschedule face neighbors after one approximate TOADS thinning removal.

This is the neighbor loop at the end of VandyCRUISE
``MultipleToadSegmentation.propagateThinning``. The selected implementation
uses the popped priority minus the reciprocal sparse-membership ratio, and
retains the source's literal 26-neighbor count condition (which tests the
outer face offset, so all 26 voxels are counted).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_heap import ToadsHeapItem, ToadsPriorityHeap4D


@dataclass(frozen=True)
class ToadsApproximateNeighborSchedule:
    """Queue additions and updated availability after one removed voxel."""

    available: NDArray[np.bool_]
    inserted: tuple[ToadsHeapItem, ...]


def schedule_approximate_thinning_neighbors(
    segmentation: ArrayLike,
    sparse_memberships: ArrayLike,
    best_class_indices: ArrayLike,
    available: ArrayLike,
    queue: ToadsPriorityHeap4D,
    removed_voxel: tuple[int, int, int],
    removed_label: int,
    popped_priority: float,
    class_volume: float,
    minimum_class_volume: float,
    limit: float,
    *,
    zero: float = 1e-10,
    copy_state: bool = True,
) -> ToadsApproximateNeighborSchedule:
    """Add the selected source's eligible six-face queue entries.

    Inputs describe the state *after* the popped voxel was removed from its
    skeleton and after its sparse-membership volume decrement. The heap is
    mutated in place, matching the selected Java tree behavior.
    """
    labels = np.asarray(segmentation)
    mems = np.asarray(sparse_memberships, dtype=np.float32)
    best = np.asarray(best_class_indices)
    state = np.asarray(available)
    if labels.ndim != 3 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be a 3-D integer grid")
    if mems.ndim != 4 or mems.shape[:3] != labels.shape or mems.shape[-1] < 2:
        raise ValueError("sparse_memberships must match segmentation and contain multiple slots")
    if best.shape != mems.shape or not np.issubdtype(best.dtype, np.integer):
        raise ValueError("best_class_indices must match sparse_memberships")
    if state.ndim != 4 or state.shape[:3] != labels.shape or state.dtype != np.bool_:
        raise ValueError("available must be a boolean grid matching segmentation and classes")
    classes = state.shape[-1]
    if classes < 2 or np.any(best < 0) or np.any(best >= classes):
        raise ValueError("sparse class indices must contain valid classes")
    if state.shape != labels.shape + (classes,):
        raise ValueError("available must be a boolean grid matching segmentation and classes")
    if not np.isfinite(mems).all():
        raise ValueError("sparse memberships must be finite")
    if not isinstance(queue, ToadsPriorityHeap4D) or queue.mode != "max":
        raise ValueError("queue must be a TOADS max heap")
    if not isinstance(removed_label, (int, np.integer)) or not 0 <= removed_label < classes:
        raise ValueError("removed_label must be a valid class index")
    x, y, z = (int(v) for v in removed_voxel)
    if not (1 <= x < labels.shape[0] - 1 and 1 <= y < labels.shape[1] - 1 and 1 <= z < labels.shape[2] - 1):
        raise ValueError("removed_voxel must be interior")
    value, volume, minimum, limit32, zero32 = map(np.float32, (popped_priority, class_volume, minimum_class_volume, limit, zero))
    if not all(np.isfinite(v) for v in (value, volume, minimum, limit32, zero32)) or zero32 <= 0 or limit32 < 0:
        raise ValueError("priority, volumes, limit, and zero must be finite; zero positive and limit nonnegative")

    updated_available = state.copy() if copy_state else state
    inserted: list[ToadsHeapItem] = []
    # Source nested i/j/l loop order, restricted to the six face neighbors.
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for l in (-1, 0, 1):
                if i * i + j * j + l * l != 1:
                    continue
                xi, yj, zl = x + i, y + j, z + l
                if not updated_available[xi, yj, zl, removed_label]:
                    continue
                # The Java inner count uses the nonzero outer (i,j,l) test.
                neighbor_count = int(np.count_nonzero(
                    labels[xi - 1:xi + 2, yj - 1:yj + 2, zl - 1:zl + 2] == removed_label
                ))
                if neighbor_count == 1:
                    order = value
                else:
                    top_class = int(best[xi, yj, zl, 0])
                    if removed_label == top_class:
                        ratio = np.float32(1.0)
                    else:
                        current = np.float32(0.0)
                        for slot in range(1, mems.shape[-1]):
                            if int(best[xi, yj, zl, slot]) == removed_label:
                                current = mems[xi, yj, zl, slot]
                        ratio = np.float32(
                            np.float32(mems[xi, yj, zl, 0] + zero32)
                            / np.float32(current + zero32)
                        )
                    with np.errstate(divide="ignore", invalid="ignore"):
                        order = np.float32(value - np.float32(np.float32(1.0) / ratio))
                if volume > minimum and order > np.float32(-limit32):
                    updated_available[xi, yj, zl, removed_label] = False
                    queue.add(order, xi, yj, zl, int(removed_label))
                    inserted.append(ToadsHeapItem(order, xi, yj, zl, int(removed_label)))
    return ToadsApproximateNeighborSchedule(updated_available, tuple(inserted))
