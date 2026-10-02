"""Sparse-membership startup for the selected approximate TOADS branch.

Implements the ``useApprox`` arm of
``MultipleToadSegmentation.initThinningLabels``. This state constructor is
deliberately separate from the exact-membership initializer because the
source stores the top-k class indices and their memberships in parallel
arrays. It does not implement a full approximate TOADS iteration.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_heap import ToadsPriorityHeap4D


@dataclass(frozen=True)
class ToadsApproximateThinningStartup:
    """Skeleton, candidate availability, sparse class volumes, and queue."""

    skeleton: NDArray[np.int16]
    available: NDArray[np.bool_]
    class_volumes: NDArray[np.float32]
    boundary_counts: NDArray[np.int32]
    queue: ToadsPriorityHeap4D


def initialize_toads_approximate_thinning(
    segmentation: ArrayLike,
    sparse_memberships: ArrayLike,
    best_class_indices: ArrayLike,
    *,
    class_count: int,
    empty_label: int = -1,
    initial_skeleton: ArrayLike | None = None,
    initial_available: ArrayLike | None = None,
    zero: float = 1e-30,
) -> ToadsApproximateThinningStartup:
    """Initialize selected ``initThinningLabels`` sparse/top-k state.

    ``sparse_memberships[..., a]`` belongs to class
    ``best_class_indices[..., a]``. The source scans interior voxels in
    x/y/z order, copies segmentation into the skeleton, accumulates only
    represented memberships for each current label, and queues boundary
    voxels by the negative reciprocal of the source top-k membership ratio.
    The 26-neighbor endpoint count and six-face boundary check preserve the
    source's distinct neighborhood semantics.
    """
    labels = np.asarray(segmentation)
    mems = np.asarray(sparse_memberships, dtype=np.float32)
    best = np.asarray(best_class_indices)
    if labels.ndim != 3 or any(length < 3 for length in labels.shape) or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer 3-D grid with each dimension at least 3")
    if mems.ndim != 4 or mems.shape[:3] != labels.shape or mems.shape[-1] < 2:
        raise ValueError("sparse_memberships must be a matching 4-D grid with at least two slots")
    if best.shape != mems.shape or not np.issubdtype(best.dtype, np.integer):
        raise ValueError("best_class_indices must be an integer grid matching sparse_memberships")
    if isinstance(class_count, bool) or not isinstance(class_count, (int, np.integer)) or class_count < 2:
        raise ValueError("class_count must be an integer of at least two")
    classes, approx = int(class_count), mems.shape[-1]
    if np.any(labels < 0) or np.any(labels >= classes) or np.any(best < 0) or np.any(best >= classes):
        raise ValueError("segmentation and sparse class indices must be nonnegative class indices")
    if np.any(labels >= classes):
        raise ValueError("segmentation labels must be valid class indices")
    if np.any(~np.isfinite(mems)) or np.any(mems < 0):
        raise ValueError("sparse memberships must be finite and nonnegative")
    if np.any(np.sort(best, axis=-1)[..., 1:] == np.sort(best, axis=-1)[..., :-1]):
        raise ValueError("sparse class indices must be unique per voxel")
    zero32 = np.float32(zero)
    if not np.isfinite(zero32) or zero32 <= np.float32(0):
        raise ValueError("zero must be finite and positive")

    if initial_skeleton is None:
        skeleton = np.full(labels.shape, empty_label, dtype=np.int16)
    else:
        source = np.asarray(initial_skeleton)
        if source.shape != labels.shape or not np.issubdtype(source.dtype, np.integer):
            raise ValueError("initial_skeleton must be an integer grid matching segmentation")
        skeleton = source.astype(np.int16, copy=True)
    if initial_available is None:
        available = np.zeros(labels.shape + (classes,), dtype=np.bool_)
    else:
        source = np.asarray(initial_available)
        if source.shape != labels.shape + (classes,) or source.dtype != np.bool_:
            raise ValueError("initial_available must be a boolean grid matching segmentation and classes")
        available = source.copy()

    nx, ny, nz = labels.shape
    volumes = np.zeros(classes, dtype=np.float32)
    boundaries = np.zeros(classes, dtype=np.int32)
    queue = ToadsPriorityHeap4D(nx * ny + ny * nz + nz * nx, "max")

    # This reproduces the source's first interior pass and float32 addition order.
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                current = int(labels[x, y, z])
                skeleton[x, y, z] = np.int16(current)
                for slot in range(approx):
                    candidate = int(best[x, y, z, slot])
                    if current == candidate:
                        volumes[current] = np.float32(volumes[current] + mems[x, y, z, slot])
                for candidate in range(classes):
                    available[x, y, z, candidate] = candidate == current

    # Source boundary selection and endpoint scoring, including strict comparisons.
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                current = int(labels[x, y, z])
                boundary = False
                neighbor_count = 0
                ni = nj = nl = 0
                for i in (-1, 0, 1):
                    for j in (-1, 0, 1):
                        for k in (-1, 0, 1):
                            d2 = i * i + j * j + k * k
                            neighbor = int(labels[x + i, y + j, z + k])
                            if d2 == 1 and neighbor != current:
                                boundary = True
                            elif d2 > 0:
                                neighbor_count += 1
                                ni, nj, nl = x + i, y + j, z + k
                if not boundary:
                    continue
                available[x, y, z, current] = False
                endpoint = False
                if neighbor_count == 1:
                    neighbor_count = 0
                    for i in (-1, 0, 1):
                        for j in (-1, 0, 1):
                            for k in (-1, 0, 1):
                                if i * i + j * j + k * k > 0 and labels[ni + i, nj + j, nl + k] == current:
                                    neighbor_count += 1
                    endpoint = neighbor_count <= 2
                if endpoint:
                    priority = np.float32(-0.0)
                else:
                    top_label = int(best[x, y, z, 0])
                    if current == top_label:
                        ratio = np.float32(1.0)
                    else:
                        current_value = np.float32(0.0)
                        for slot in range(1, approx):
                            if int(best[x, y, z, slot]) == current:
                                current_value = mems[x, y, z, slot]
                        ratio = np.float32(
                            np.float32(mems[x, y, z, 0] + zero32)
                            / np.float32(current_value + zero32)
                        )
                    priority = np.float32(np.float32(-1.0) / ratio)
                queue.add(priority, x, y, z, current)
                boundaries[current] += 1

    return ToadsApproximateThinningStartup(skeleton, available, volumes, boundaries, queue)
