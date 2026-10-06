"""Source-ordered sparse-membership TOADS thinning queue pass.

This composes the already isolated approximate constructor initializer,
selected per-item relation/topology update, and Java-order neighbor scheduler.
It is a queue-stage API; it does not calculate approximate memberships or run
the surrounding centroid, registration, growing, and convergence workflow.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_approximate_neighbor_schedule import schedule_approximate_thinning_neighbors
from .toads_approximate_thinning_startup import (
    ToadsApproximateThinningStartup,
    initialize_toads_approximate_thinning,
)
from .toads_approximate_thinning_step import (
    ToadsApproximateThinningStep,
    process_next_approximate_thinning_item,
)
from .toads_heap import ToadsHeapItem, ToadsPriorityHeap4D


@dataclass(frozen=True)
class ToadsApproximateThinningPass:
    """Terminal state and compact optional event history for one queue pass."""

    segmentation: NDArray[np.int16]
    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    class_volumes: NDArray[np.float32]
    minimum_class_volumes: NDArray[np.float32]
    boundary_counts: NDArray[np.int32]
    queue: ToadsPriorityHeap4D
    processed_count: int
    stop_reason: str
    processed_items: tuple[ToadsHeapItem, ...]


def run_approximate_thinning_pass(
    segmentation: ArrayLike,
    sparse_memberships: ArrayLike,
    best_class_indices: ArrayLike,
    relation_matrix: ArrayLike,
    transformed_shape_values: ArrayLike,
    *,
    limit: float,
    volume_ratio: float = 0.0,
    minimum_volume_floor: float = 100.0,
    explicit_minimum_class_volumes: ArrayLike | None = None,
    class_count: int | None = None,
    empty_label: int = -1,
    record_processed_items: bool = False,
) -> ToadsApproximateThinningPass:
    """Initialize and drain one sparse-membership ``propagateThinning`` pass.

    The class minimum defaults to the selected Java rule
    ``max(minimum_volume_floor, volume_ratio * initial_class_volume)``. An
    explicit minimum vector supports small synthetic references that need to
    exercise queue transitions below the production 100-voxel floor.
    """
    labels_in = np.asarray(segmentation)
    memberships = np.asarray(sparse_memberships, dtype=np.float32)
    best = np.asarray(best_class_indices)
    relations = np.asarray(relation_matrix)
    shapes = np.asarray(transformed_shape_values, dtype=np.float32)
    if labels_in.ndim != 3 or not np.issubdtype(labels_in.dtype, np.integer):
        raise ValueError("segmentation must be a 3-D integer grid")
    if memberships.ndim != 4 or memberships.shape[:3] != labels_in.shape:
        raise ValueError("sparse_memberships must match segmentation")
    classes = int(class_count) if class_count is not None else int(np.max(best)) + 1
    if classes < 2:
        raise ValueError("class_count must be at least two")
    if relations.shape != (classes, classes) or relations.dtype != np.bool_:
        raise ValueError("relation_matrix must be a boolean square matrix matching classes")
    if shapes.shape != labels_in.shape + (classes,) or not np.isfinite(shapes).all():
        raise ValueError("transformed_shape_values must be finite on the segmentation grid")
    limit32, ratio32, floor32 = map(np.float32, (limit, volume_ratio, minimum_volume_floor))
    if not all(np.isfinite(v) for v in (limit32, ratio32, floor32)) or limit32 < 0 or ratio32 < 0 or floor32 < 0:
        raise ValueError("limit, volume_ratio, and minimum_volume_floor must be finite and nonnegative")
    labels = labels_in.astype(np.int16, copy=True)
    startup = initialize_toads_approximate_thinning(
        labels,
        memberships,
        best,
        class_count=classes,
        empty_label=empty_label,
    )
    if explicit_minimum_class_volumes is None:
        minimums = np.maximum(
            floor32,
            np.asarray(ratio32 * startup.class_volumes, dtype=np.float32),
        ).astype(np.float32)
    else:
        minimums = np.asarray(explicit_minimum_class_volumes, dtype=np.float32)
        if minimums.shape != (classes,) or not np.isfinite(minimums).all() or np.any(minimums < 0):
            raise ValueError("explicit_minimum_class_volumes must be a finite nonnegative class vector")
        minimums = minimums.copy()

    skeleton = startup.skeleton
    available = startup.available
    volumes = startup.class_volumes
    queue = startup.queue
    processed: list[ToadsHeapItem] = []
    processed_count = 0
    previous_priority = np.float32(0.0)
    while queue and previous_priority > np.float32(-limit32):
        current_startup = ToadsApproximateThinningStartup(
            skeleton, available, volumes, startup.boundary_counts, queue
        )
        event: ToadsApproximateThinningStep = process_next_approximate_thinning_item(
            labels,
            memberships,
            best,
            current_startup,
            relations,
            shapes,
            minimums,
            empty_label=empty_label,
            copy_state=False,
            validate_full_arrays=False,
        )
        labels = event.segmentation
        skeleton = event.skeleton
        available = event.available
        volumes = event.class_volumes
        queue = event.remaining_queue
        previous_priority = event.popped_item.value
        if record_processed_items:
            processed.append(event.popped_item)
        processed_count += 1
        if event.transition is not None and event.transition.eligible:
            scheduled = schedule_approximate_thinning_neighbors(
                labels,
                memberships,
                best,
                available,
                queue,
                (event.popped_item.x, event.popped_item.y, event.popped_item.z),
                event.popped_item.label,
                event.popped_item.value,
                volumes[event.popped_item.label],
                minimums[event.popped_item.label],
                limit32,
                copy_state=False,
                validate_full_arrays=False,
            )
            available = scheduled.available

    stop_reason = "priority_limit" if queue and previous_priority <= np.float32(-limit32) else "queue_empty"
    return ToadsApproximateThinningPass(
        labels,
        skeleton,
        available,
        volumes,
        minimums,
        startup.boundary_counts,
        queue,
        processed_count,
        stop_reason,
        tuple(processed),
    )
