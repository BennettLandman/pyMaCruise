"""One sparse-membership queue event after approximate TOADS initialization.

This bridge handles one popped ``propagateThinning`` item. The shared local
transition owns relation-aware target selection and the selected 18/6
homeomorphism gate; this adapter preserves the source's ``best[..., 0]``
target and sparse class-volume decrement. Neighbor rescheduling and a full
outer queue pass remain separate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_approximate_thinning_startup import ToadsApproximateThinningStartup
from .toads_heap import ToadsHeapItem, ToadsPriorityHeap4D
from .toads_thinning import ToadsThinningTransition, toads_thinning_voxel_transition


@dataclass(frozen=True)
class ToadsApproximateThinningStep:
    """Image state after processing one item from the approximate queue."""

    segmentation: NDArray[np.int16]
    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    class_volumes: NDArray[np.float32]
    popped_item: ToadsHeapItem
    transition: ToadsThinningTransition | None
    remaining_queue: ToadsPriorityHeap4D


def process_next_approximate_thinning_item(
    segmentation: ArrayLike,
    sparse_memberships: ArrayLike,
    best_class_indices: ArrayLike,
    startup: ToadsApproximateThinningStartup,
    relation_matrix: ArrayLike,
    transformed_shape_values: ArrayLike,
    minimum_class_volumes: ArrayLike,
    *,
    empty_label: int = -1,
    copy_state: bool = True,
    validate_full_arrays: bool = True,
) -> ToadsApproximateThinningStep:
    """Pop and apply one approximate-path ``propagateThinning`` event.

    The Java approximate branch takes its proposed tissue from sparse slot 0,
    irrespective of which sparse membership value is largest. The local
    transition shared with the exact path expects a dense class vector; this
    bridge maps represented slots to their class labels and raises slot 0's
    score only enough to preserve the selected Java target. Since that score
    is only a target-selection device, class volume is then updated using the
    original sparse membership for the popped class, exactly as Java does.

    This function deliberately stops after the selected voxel's label,
    skeleton, availability, and class-volume update. The source's six-face
    neighbor rescheduling is not included in this one-event adapter.
    """
    labels_in = np.asarray(segmentation)
    mems = np.asarray(sparse_memberships, dtype=np.float32)
    best = np.asarray(best_class_indices)
    if labels_in.ndim != 3 or not np.issubdtype(labels_in.dtype, np.integer):
        raise ValueError("segmentation must be a 3-D integer grid")
    if mems.ndim != 4 or mems.shape[:3] != labels_in.shape or mems.shape[-1] < 2:
        raise ValueError("sparse_memberships must match segmentation and contain multiple slots")
    if best.shape != mems.shape or not np.issubdtype(best.dtype, np.integer):
        raise ValueError("best_class_indices must match sparse_memberships")
    classes = startup.class_volumes.size
    if np.any(best < 0) or np.any(best >= classes) or np.any(labels_in < 0) or np.any(labels_in >= classes):
        raise ValueError("segmentation and sparse class labels must be valid startup classes")
    if startup.skeleton.shape != labels_in.shape or startup.available.shape != labels_in.shape + (classes,):
        raise ValueError("startup state must match segmentation and class dimensions")
    relations = np.asarray(relation_matrix)
    shapes = np.asarray(transformed_shape_values, dtype=np.float32)
    minima = np.asarray(minimum_class_volumes, dtype=np.float32)
    if relations.shape != (classes, classes) or relations.dtype != np.bool_:
        raise ValueError("relation_matrix must be a boolean square matrix matching classes")
    if shapes.shape != labels_in.shape + (classes,) or not np.isfinite(shapes).all():
        raise ValueError("transformed_shape_values must be finite on the segmentation grid")
    if minima.shape != (classes,) or not np.isfinite(minima).all():
        raise ValueError("minimum_class_volumes must be finite and match classes")
    if validate_full_arrays and not np.isfinite(mems).all():
        raise ValueError("sparse memberships must be finite")
    if not startup.queue:
        raise ValueError("cannot process an empty approximate thinning queue")

    updated_segmentation = labels_in.astype(np.int16, copy=copy_state)
    updated_skeleton = startup.skeleton.copy() if copy_state else startup.skeleton
    updated_available = startup.available.copy() if copy_state else startup.available
    updated_volumes = startup.class_volumes.copy() if copy_state else startup.class_volumes
    popped = startup.queue.pop_first()
    x, y, z, processing_label = popped.x, popped.y, popped.z, popped.label
    transition: ToadsThinningTransition | None = None
    if (
        1 <= x < labels_in.shape[0] - 1
        and 1 <= y < labels_in.shape[1] - 1
        and 1 <= z < labels_in.shape[2] - 1
        and updated_skeleton[x, y, z] == processing_label
    ):
        dense = np.zeros(classes, dtype=np.float32)
        for slot in range(mems.shape[-1]):
            dense[int(best[x, y, z, slot])] = mems[x, y, z, slot]
        top_class = int(best[x, y, z, 0])
        if dense[top_class] <= np.max(dense):
            dense[top_class] = np.nextafter(np.max(dense), np.float32(np.inf), dtype=np.float32)
        neighborhood = updated_segmentation[x - 1:x + 2, y - 1:y + 2, z - 1:z + 2]
        skeleton_neighborhood = updated_skeleton[x - 1:x + 2, y - 1:y + 2, z - 1:z + 2]
        shape_neighborhood = shapes[x - 1:x + 2, y - 1:y + 2, z - 1:z + 2]
        transition = toads_thinning_voxel_transition(
            neighborhood,
            skeleton_neighborhood,
            dense,
            relations,
            shape_neighborhood,
            updated_volumes,
            minima,
            available_at_center=updated_available[x, y, z],
            processing_label=processing_label,
            empty_label=empty_label,
        )
        if transition.eligible:
            updated_segmentation[x, y, z] = np.int16(transition.segmentation_label)
            updated_skeleton[x, y, z] = np.int16(transition.skeleton_label)
            updated_available[x, y, z] = transition.available_at_center
            sparse_membership = np.float32(0.0)
            for slot in range(mems.shape[-1]):
                if int(best[x, y, z, slot]) == processing_label:
                    sparse_membership = np.float32(sparse_membership + mems[x, y, z, slot])
            updated_volumes[processing_label] = np.float32(
                updated_volumes[processing_label] - sparse_membership
            )

    return ToadsApproximateThinningStep(
        updated_segmentation,
        updated_skeleton,
        updated_available,
        updated_volumes,
        popped,
        transition,
        startup.queue,
    )
