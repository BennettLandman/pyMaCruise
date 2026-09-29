"""Topology-preserving thinning from the supplied VandyCRUISE ThinACE source.

This is a separate source-specific cleanup stage after ACE/MaACE has found
candidate separating voxels. Tosun et al. (2004) and Huo et al. (2016) describe
thin-skeleton output as part of ACE/MaACE, but the priority queue, simple-point
tests, connectivity conventions, and endpoint rules here are traced to the
VandyCRUISE ``ThinACE`` implementation rather than inferred from those papers.
See ``docs/papers/ace-maace.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib
import numpy as np
from scipy import ndimage

from pymacruise.ace import _as_ras_3d, _make_output, _signed_distance_fmm

_OFFSETS_6 = tuple(
    (dx, dy, dz)
    for dx, dy, dz in (
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (-1, 0, 0),
        (0, -1, 0),
        (0, 0, -1),
    )
)
_OFFSETS_SURFACE_EXPANSION = (
    (0, 1, 0),
    (1, 0, 0),
    (0, 0, 1),
    (0, -1, 0),
    (-1, 0, 0),
    (0, 0, -1),
)
_OFFSETS_18 = tuple(
    (dx, dy, dz)
    for dx in (-1, 0, 1)
    for dy in (-1, 0, 1)
    for dz in (-1, 0, 1)
    if 0 < abs(dx) + abs(dy) + abs(dz) <= 2
)
_OFFSETS_26 = tuple(
    (dx, dy, dz)
    for dx in (-1, 0, 1)
    for dy in (-1, 0, 1)
    for dz in (-1, 0, 1)
    if (dx, dy, dz) != (0, 0, 0)
)
_FOREGROUND18 = ndimage.generate_binary_structure(3, 2)
_FOREGROUND18_DOMAIN = np.zeros((3, 3, 3), dtype=bool)
for _dx, _dy, _dz in _OFFSETS_18:
    _FOREGROUND18_DOMAIN[_dx + 1, _dy + 1, _dz + 1] = True
_BACKGROUND6_DOMAIN = np.zeros((3, 3, 3), dtype=bool)
for _dx, _dy, _dz in _OFFSETS_6:
    _BACKGROUND6_DOMAIN[_dx + 1, _dy + 1, _dz + 1] = True


@dataclass(frozen=True)
class ThinACEImageResult:
    """RAS-oriented NIfTI result of the source ThinACE operation."""

    skeleton: nib.Nifti1Image


def thin_ace(skeleton: np.ndarray) -> np.ndarray:
    """Thin a positive ACE skeleton with distance-priority simple-point peeling.

    This follows the source's iteration structure, 18-connected foreground /
    6-connected background tests, endpoint retention, and zero image border.
    Output foreground is encoded as 254, matching `ThinACE.solve`.

    The papers motivate thinning the ACE shock set to a one-voxel skeleton;
    these precise topology and ordering semantics are the VandyCRUISE source
    implementation and should not be attributed to the paper-level method.
    """
    original = np.asarray(skeleton)
    if original.ndim != 3 or any(size < 3 for size in original.shape):
        raise ValueError("ThinACE requires a 3D array with dimensions >= 3")
    if not np.isfinite(original).all():
        raise ValueError("ThinACE input must be finite")

    object_mask = original > 0
    initial_phi = np.where(object_mask, 1.0, -1.0).astype(np.float32)
    object_mask = object_mask.copy()
    object_mask[[0, -1], :, :] = False
    object_mask[:, [0, -1], :] = False
    object_mask[:, :, [0, -1]] = False
    distance = _signed_distance_fmm(
        initial_phi,
        max_distance=4.0,
        round_updates_to_float32=True,
    )

    labels = np.zeros(object_mask.shape, dtype=np.uint8)
    list1 = [
        (i, j, k)
        for k in range(1, object_mask.shape[2] - 1)
        for j in range(1, object_mask.shape[1] - 1)
        for i in range(1, object_mask.shape[0] - 1)
        if object_mask[i, j, k]
        and any(not object_mask[i + dx, j + dy, k + dz] for dx, dy, dz in _OFFSETS_6)
    ]
    if list1:
        labels[tuple(np.asarray(list1, dtype=np.intp).T)] = 1

    while list1:
        list2: list[tuple[int, int, int]] = []
        heap: list[tuple[np.float32, tuple[int, int, int]]] = []
        found_simple = False

        for point in list1:
            if not object_mask[point]:
                continue
            foreground_components, background_components = _simple_point_counts(
                object_mask, point
            )
            if foreground_components == 1 and background_components == 2:
                labels[point] = 4
            elif foreground_components == 1 and background_components == 1:
                i, j, k = point
                _source_heap_add(heap, np.float32(distance[point]), point)
                found_simple = True
            else:
                labels[point] = 2
                list2.append(point)

        while heap:
            _, point = _source_heap_remove(heap)
            i, j, k = point
            labels[point] = 0
            object_neighbors = _neighbor_values(object_mask, point, _OFFSETS_26)
            count26 = sum(object_neighbors)
            count18 = sum(_neighbor_values(object_mask, point, _OFFSETS_18))
            foreground_components, background_components = _simple_point_counts(
                object_mask, point
            )
            if foreground_components != 1 or background_components != 1:
                labels[point] = 2
                list2.append(point)
                continue

            if count18 == 1:
                count26 = 1
            if _is_source_endpoint(object_mask, point, count26):
                labels[point] = 4
                continue

            object_mask[point] = False
            for dx, dy, dz in _OFFSETS_SURFACE_EXPANSION:
                neighbor = (i + dx, j + dy, k + dz)
                if object_mask[neighbor] and labels[neighbor] == 0:
                    labels[neighbor] = 2
                    list2.append(neighbor)

        if not found_simple:
            break
        list1 = list2

    return (object_mask.astype(np.uint8) * 254).astype(np.int32)


def thin_ace_image(skeleton: nib.spatialimages.SpatialImage) -> ThinACEImageResult:
    """Run ThinACE on a NIfTI image, warning and reorienting to RAS if needed."""
    image = _as_ras_3d(skeleton, "ACE skeleton")
    result = thin_ace(np.asarray(image.dataobj))
    return ThinACEImageResult(_make_output(result, image))


def _neighbor_values(
    object_mask: np.ndarray,
    point: tuple[int, int, int],
    offsets: tuple[tuple[int, int, int], ...],
) -> list[bool]:
    i, j, k = point
    return [bool(object_mask[i + dx, j + dy, k + dz]) for dx, dy, dz in offsets]


def _simple_point_counts(
    object_mask: np.ndarray,
    point: tuple[int, int, int],
) -> tuple[int, int]:
    neighbors = object_mask[
        point[0] - 1 : point[0] + 2,
        point[1] - 1 : point[1] + 2,
        point[2] - 1 : point[2] + 2,
    ]
    # ThinACE does not treat every non-object location as background. It
    # propagates the negative label outward from 6-neighbor background to
    # 18-neighbors, then to corners; untouched locations remain neutral.
    cube = np.zeros((3, 3, 3), dtype=np.int8)
    for dx, dy, dz in _OFFSETS_6:
        cube[dx + 1, dy + 1, dz + 1] = 1 if neighbors[dx + 1, dy + 1, dz + 1] else -1
    for dx, dy, dz in _OFFSETS_18:
        index = (dx + 1, dy + 1, dz + 1)
        if neighbors[index]:
            cube[index] = 1
        elif dz == 0:
            if cube[1, dy + 1, 1] == -1 or cube[dx + 1, 1, 1] == -1:
                cube[index] = -1
        elif dy == 0:
            if cube[dx + 1, 1, 1] == -1 or cube[1, 1, dz + 1] == -1:
                cube[index] = -1
        elif cube[1, 1, dz + 1] == -1 or cube[1, dy + 1, 1] == -1:
            cube[index] = -1
    for dx in (-1, 1):
        for dy in (-1, 1):
            for dz in (-1, 1):
                index = (dx + 1, dy + 1, dz + 1)
                if neighbors[index]:
                    cube[index] = 1
                elif (
                    cube[dx + 1, dy + 1, 1] == -1
                    or cube[dx + 1, 1, dz + 1] == -1
                    or cube[1, dy + 1, dz + 1] == -1
                ):
                    cube[index] = -1

    # The legacy component routines enumerate only their labeled neighborhood
    # domain (18 neighbors for foreground; six for background); disconnected
    # values outside that domain are intentionally ignored.
    foreground = cube == 1
    background = cube == -1
    # The source BFS stops after discovering a second component and returns 2
    # as a sentinel, so preserve that cap for its special surface predicate.
    return (
        _components_touching_domain(foreground, _FOREGROUND18, _FOREGROUND18_DOMAIN),
        _components_touching_domain(background, ndimage.generate_binary_structure(3, 1), _BACKGROUND6_DOMAIN),
    )


def _components_touching_domain(
    mask: np.ndarray,
    structure: np.ndarray,
    seed_domain: np.ndarray,
) -> int:
    labels, _ = ndimage.label(mask, structure=structure)
    components = np.unique(labels[seed_domain])
    count = int(np.count_nonzero(components))
    # Source connectivity search begins only at the requested neighbor domain
    # but traverses the full 3x3x3 cube; it returns as soon as a second
    # connected component is found.
    return min(count, 2)


def _source_heap_add(
    heap: list[tuple[np.float32, tuple[int, int, int]]],
    distance: np.float32,
    point: tuple[int, int, int],
) -> None:
    """Reproduce JIST BinaryMinHeap's strict-comparison tie behavior."""
    hole = len(heap)
    heap.append((distance, point))
    while hole > 0:
        parent = (hole - 1) // 2
        if not distance < heap[parent][0]:
            break
        heap[hole] = heap[parent]
        hole = parent
    heap[hole] = (distance, point)


def _source_heap_remove(
    heap: list[tuple[np.float32, tuple[int, int, int]]],
) -> tuple[np.float32, tuple[int, int, int]]:
    """Remove the min using the child/tie rules in JIST BinaryMinHeap."""
    minimum = heap[0]
    last = heap.pop()
    if heap:
        hole = 0
        while 2 * hole + 1 < len(heap):
            child = 2 * hole + 1
            right = child + 1
            if right < len(heap) and heap[right][0] < heap[child][0]:
                child = right
            if heap[child][0] < last[0]:
                heap[hole] = heap[child]
                hole = child
            else:
                break
        heap[hole] = last
    return minimum


def _is_source_endpoint(
    object_mask: np.ndarray,
    point: tuple[int, int, int],
    count26: int,
) -> bool:
    if count26 in (2, 3):
        return True
    if not 4 <= count26 <= 7:
        return False

    i, j, k = point
    for sx in (-1, 1):
        for sy in (-1, 1):
            for sz in (-1, 1):
                octant = (
                    (sx, 0, 0),
                    (0, sy, 0),
                    (0, 0, sz),
                    (sx, sy, 0),
                    (sx, 0, sz),
                    (0, sy, sz),
                    (sx, sy, sz),
                )
                occupied = sum(
                    object_mask[i + dx, j + dy, k + dz]
                    for dx, dy, dz in octant
                )
                if occupied > 3:
                    return False
    return True
