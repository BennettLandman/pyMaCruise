"""Deterministic constructor state for the selected LesionTOADS engine.

The selected VandyCRUISE ``LesionToadSegmentation`` constructor copies its
atlas classification, maps template labels to class indices, accumulates
class volumes and channel centroids, and creates initial zero/one model
buffers before segmentation evolution begins. This module isolates those
constructor effects. It does not estimate memberships or evolve topology.

The numeric loop order follows the Java ``x, y, z, class`` traversal. Arrays
enter this helper in the package's RAS voxel-axis convention; image loading
and orientation warnings remain the caller's responsibility.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class LesionToadsSegmentationState:
    """Initial arrays and statistics created by the Java constructor."""

    segmentation: NDArray[np.int8]
    skeleton: NDArray[np.int8]
    relation_map: NDArray[np.int16]
    memberships: NDArray[np.float32]
    available: NDArray[np.bool_]
    class_volume: NDArray[np.float32]
    image_centroid: NDArray[np.float32]
    intensity_stddev: NDArray[np.float32]
    ordering: NDArray[np.float32]
    lesion_classification: NDArray[np.int16]
    lesion_classification_membership: NDArray[np.int16]
    lesion_membership: NDArray[np.float32]
    lesion_relation: NDArray[np.float32]
    outlier_membership: NDArray[np.float32]
    black_hole_membership: NDArray[np.float32]
    ventricle_distance: NDArray[np.int8]
    gray_matter_distance: NDArray[np.int8]
    brainstem_distance: NDArray[np.int8]
    csf_distance: NDArray[np.int8]
    ventricle_mask: NDArray[np.bool_]
    gray_matter_mask: NDArray[np.bool_]
    brainstem_mask: NDArray[np.bool_]
    csf_mask: NDArray[np.bool_]
    covariance: NDArray[np.float32]
    lesion_covariance: NDArray[np.float32]
    class_similarity: NDArray[np.float32]
    lesion_similarity: NDArray[np.float32]
    black_hole_similarity: NDArray[np.float32]
    prior_centroid: NDArray[np.float32]
    lesion_centroid: NDArray[np.float32]
    black_hole_centroid: NDArray[np.float32]
    lesion_stddev: NDArray[np.float32]
    black_hole_stddev: NDArray[np.float32]
    image_high: NDArray[np.float32]
    image_scale: NDArray[np.float32]
    channel_variance: NDArray[np.float32]
    channel_weight: NDArray[np.float32]
    lesion_weight: NDArray[np.float32]
    lesion_class: int


def initialize_lesion_toads_segmentation_state(
    images: ArrayLike,
    template_labels: ArrayLike,
    class_labels: Sequence[int],
    object_types: Sequence[str],
    *,
    max_ventricle_distance: int = 9,
    max_gray_matter_distance: int = 5,
    max_brainstem_distance: int = 5,
) -> LesionToadsSegmentationState:
    """Create the selected Java constructor's initial segmentation buffers.

    ``images`` is ``(channels, x, y, z)`` and ``template_labels`` is
    ``(x, y, z)``. ``class_labels`` and ``object_types`` use atlas class
    order. Class centroids are accumulated with float32 updates in source
    traversal order. An absent class has the Java ``0/0`` NaN centroid.

    The Java constructor casts configured distance maxima to signed bytes,
    then stores the maximum plus one (also as a signed byte). This helper
    retains that behavior. Unknown template labels reject explicitly rather
    than returning the source's partially initialized, aborted object.
    """
    image_values = np.asarray(images, dtype=np.float32)
    labels = np.asarray(template_labels)
    atlas_labels = np.asarray(class_labels, dtype=np.int64)
    types = tuple(str(value) for value in object_types)

    if image_values.ndim != 4 or image_values.shape[0] < 1:
        raise ValueError("images must have shape (channels, x, y, z)")
    if labels.ndim != 3 or any(size == 0 for size in labels.shape):
        raise ValueError("template_labels must be a nonempty 3D array")
    if image_values.shape[1:] != labels.shape:
        raise ValueError("image and template-label grids must match")
    if not np.isfinite(image_values).all():
        raise ValueError("images must contain only finite values")
    if not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("template labels must be integers")
    if atlas_labels.ndim != 1 or atlas_labels.size == 0:
        raise ValueError("at least one ordered class label is required")
    if atlas_labels.size > 128:
        raise ValueError("class indices must fit the Java signed-byte storage")
    if len(types) != atlas_labels.size:
        raise ValueError("one object type is required for each class label")
    if np.any(atlas_labels < -128) or np.any(atlas_labels > 127):
        raise ValueError("class labels must fit the Java signed-byte template")
    if np.any(labels < -128) or np.any(labels > 127):
        raise ValueError("template voxel labels must fit the Java signed-byte template")
    if not all(isinstance(value, str) and value for value in object_types):
        raise ValueError("object types must be nonempty strings")
    maxima = (
        max_ventricle_distance,
        max_gray_matter_distance,
        max_brainstem_distance,
    )
    if any(not isinstance(value, (int, np.integer)) for value in maxima):
        raise ValueError("distance maxima must be integers")
    if any(value < -32768 or value > 32767 for value in maxima):
        raise ValueError("distance maxima must fit the Java short constructor inputs")

    class_count = int(atlas_labels.size)
    channel_count = int(image_values.shape[0])
    grid = tuple(int(size) for size in labels.shape)
    segmentation = np.full(grid, -1, dtype=np.int8)
    skeleton = np.full(grid, -1, dtype=np.int8)
    relation_map = np.full(grid, -1, dtype=np.int16)
    memberships = np.zeros(grid + (class_count,), dtype=np.float32)
    available = np.zeros(grid + (class_count,), dtype=np.bool_)
    class_volume = np.zeros(class_count, dtype=np.float32)
    image_centroid = np.zeros((channel_count, class_count), dtype=np.float32)

    # Match Java's nested x/y/z/class loops and single-precision += order.
    for x in range(grid[0]):
        for y in range(grid[1]):
            for z in range(grid[2]):
                matched = False
                voxel_label = int(labels[x, y, z])
                for class_index, atlas_label in enumerate(atlas_labels):
                    if voxel_label == int(atlas_label):
                        matched = True
                        class_volume[class_index] = np.float32(
                            class_volume[class_index] + np.float32(1.0)
                        )
                        segmentation[x, y, z] = np.int8(class_index)
                        skeleton[x, y, z] = np.int8(class_index)
                        relation_map[x, y, z] = np.int16(class_index)
                        for channel in range(channel_count):
                            image_centroid[channel, class_index] = np.float32(
                                image_centroid[channel, class_index]
                                + image_values[channel, x, y, z]
                            )
                if not matched:
                    raise ValueError(
                        "template contains a label absent from class_labels: "
                        f"{voxel_label} at {(x, y, z)}"
                    )

    with np.errstate(divide="ignore", invalid="ignore"):
        for channel in range(channel_count):
            for class_index in range(class_count):
                image_centroid[channel, class_index] = np.float32(
                    image_centroid[channel, class_index]
                    / class_volume[class_index]
                )

    intensity_stddev = np.ones((channel_count, class_count), dtype=np.float32)
    ordering = np.zeros(grid, dtype=np.float32)
    lesion_classification = np.zeros(grid, dtype=np.int16)
    lesion_classification_membership = np.zeros(grid, dtype=np.int16)
    lesion_membership = np.zeros(grid, dtype=np.float32)
    lesion_relation = np.ones(grid, dtype=np.float32)
    outlier_membership = np.full(grid, np.float32(0.5), dtype=np.float32)
    black_hole_membership = np.zeros(grid, dtype=np.float32)

    def java_byte(value: int) -> np.int8:
        return np.int8(((value + 128) % 256) - 128)

    def distance_sentinel(value: int) -> np.int8:
        # The source first casts the short parameter to byte, then adds one.
        configured = int(java_byte(value))
        return java_byte(configured + 1)

    ventricle_distance = np.full(
        grid, distance_sentinel(max_ventricle_distance), dtype=np.int8
    )
    gray_matter_distance = np.full(
        grid, distance_sentinel(max_gray_matter_distance), dtype=np.int8
    )
    brainstem_distance = np.full(
        grid, distance_sentinel(max_brainstem_distance), dtype=np.int8
    )
    csf_distance = np.zeros(grid, dtype=np.int8)
    ventricle_mask = np.zeros(grid, dtype=np.bool_)
    gray_matter_mask = np.zeros(grid, dtype=np.bool_)
    brainstem_mask = np.zeros(grid, dtype=np.bool_)
    csf_mask = np.zeros(grid, dtype=np.bool_)

    covariance = np.zeros((class_count, channel_count, channel_count), dtype=np.float32)
    for class_index in range(class_count):
        covariance[class_index] = np.eye(channel_count, dtype=np.float32)
    lesion_covariance = np.eye(channel_count, dtype=np.float32)
    class_similarity = np.ones((class_count, class_count), dtype=np.float32)
    np.fill_diagonal(class_similarity, np.float32(0.0))
    lesion_similarity = np.ones(class_count, dtype=np.float32)
    black_hole_similarity = np.ones(class_count, dtype=np.float32)

    lesion_types = {"lesion", "blackhole", "optimized_lesion"}
    has_lesion = any(value in lesion_types or value == "optimized_blackHole" for value in types)
    lesion_class = -1
    if has_lesion:
        for class_index, object_type in enumerate(types):
            if object_type in lesion_types:
                lesion_class = class_index
        if lesion_class >= 0:
            lesion_similarity[lesion_class] = np.float32(0.0)
            black_hole_similarity[lesion_class] = np.float32(0.0)

    prior_centroid = np.zeros((channel_count, class_count), dtype=np.float32)
    lesion_centroid = np.zeros(channel_count, dtype=np.float32)
    black_hole_centroid = np.zeros(channel_count, dtype=np.float32)
    lesion_stddev = np.ones(channel_count, dtype=np.float32)
    black_hole_stddev = np.ones(channel_count, dtype=np.float32)
    image_high = np.ones(channel_count, dtype=np.float32)
    image_scale = np.ones(channel_count, dtype=np.float32)
    channel_variance = np.ones(channel_count, dtype=np.float32)
    channel_weight = np.full(
        channel_count, np.float32(1.0 / channel_count), dtype=np.float32
    )
    lesion_weight = channel_weight.copy()

    return LesionToadsSegmentationState(
        segmentation=segmentation,
        skeleton=skeleton,
        relation_map=relation_map,
        memberships=memberships,
        available=available,
        class_volume=class_volume,
        image_centroid=image_centroid,
        intensity_stddev=intensity_stddev,
        ordering=ordering,
        lesion_classification=lesion_classification,
        lesion_classification_membership=lesion_classification_membership,
        lesion_membership=lesion_membership,
        lesion_relation=lesion_relation,
        outlier_membership=outlier_membership,
        black_hole_membership=black_hole_membership,
        ventricle_distance=ventricle_distance,
        gray_matter_distance=gray_matter_distance,
        brainstem_distance=brainstem_distance,
        csf_distance=csf_distance,
        ventricle_mask=ventricle_mask,
        gray_matter_mask=gray_matter_mask,
        brainstem_mask=brainstem_mask,
        csf_mask=csf_mask,
        covariance=covariance,
        lesion_covariance=lesion_covariance,
        class_similarity=class_similarity,
        lesion_similarity=lesion_similarity,
        black_hole_similarity=black_hole_similarity,
        prior_centroid=prior_centroid,
        lesion_centroid=lesion_centroid,
        black_hole_centroid=black_hole_centroid,
        lesion_stddev=lesion_stddev,
        black_hole_stddev=black_hole_stddev,
        image_high=image_high,
        image_scale=image_scale,
        channel_variance=channel_variance,
        channel_weight=channel_weight,
        lesion_weight=lesion_weight,
        lesion_class=lesion_class,
    )
