"""Initialize selected TOADS constructor arrays from transformed atlas labels.

The selected ``MultipleToadSegmentation`` constructor maps each transformed
template label to its descriptor-order class index and seeds segmentation,
skeleton, and relation maps from that index. It also initializes class
volumes and zero boundary distances. These are wrapper/constructor state
details; topology labels and all outputs are kept in the caller's RAS grid.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class ToadsSegmentationConstructorState:
    """Initial array state allocated by the selected segmentation constructor."""

    classification_labels: NDArray
    segmentation: NDArray[np.int8]
    skeleton: NDArray[np.int8]
    relation_map: NDArray[np.int16]
    available: NDArray[np.bool_]
    class_volumes: NDArray[np.float32]
    boundary_distances: NDArray[np.int8]
    call_order: tuple[str, ...]


TOADS_SEGMENTATION_CONSTRUCTOR_CALL_ORDER = (
    "segmentation.constructor.mapTemplateLabels",
    "segmentation.constructor.initializeSegmentationAndSkeleton",
    "segmentation.constructor.initializeRelationAndBoundaryState",
)


def initialize_toads_segmentation_constructor_state(
    transformed_topology: Any,
    class_labels: Any,
) -> ToadsSegmentationConstructorState:
    """Map source topology labels to source class indices and initial arrays."""
    labels = np.asarray(transformed_topology)
    class_values = np.asarray(class_labels)
    if labels.ndim != 3 or any(size <= 0 for size in labels.shape):
        raise ValueError("transformed_topology must be a nonempty 3-D label map")
    if not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("transformed_topology must contain integer labels")
    if class_values.ndim != 1 or class_values.size < 1 or not np.issubdtype(class_values.dtype, np.integer):
        raise ValueError("class_labels must be a nonempty integer vector in descriptor order")
    if class_values.size > 127:
        raise ValueError("TOADS class count must fit source signed-byte class indices")
    if np.unique(class_values).size != class_values.size:
        raise ValueError("class_labels must be unique")

    segmentation = np.full(labels.shape, -1, dtype=np.int8)
    for class_index, source_label in enumerate(class_values):
        segmentation[labels == source_label] = np.int8(class_index)
    if np.any(segmentation < 0):
        unknown = np.unique(labels[segmentation < 0]).tolist()
        raise ValueError(f"transformed_topology contains labels absent from the atlas descriptor: {unknown}")

    class_count = int(class_values.size)
    counts = np.bincount(segmentation.ravel().astype(np.intp), minlength=class_count)
    volumes = counts.astype(np.float32)
    return ToadsSegmentationConstructorState(
        classification_labels=labels.copy(),
        segmentation=segmentation,
        skeleton=segmentation.copy(),
        relation_map=segmentation.astype(np.int16),
        available=np.zeros(labels.shape + (class_count,), dtype=bool),
        class_volumes=volumes,
        boundary_distances=np.zeros(labels.shape + (class_count,), dtype=np.int8),
        call_order=TOADS_SEGMENTATION_CONSTRUCTOR_CALL_ORDER,
    )
