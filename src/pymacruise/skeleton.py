"""Topology-preserving thinning from the supplied VandyCRUISE ThinACE source."""

from __future__ import annotations

import heapq
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
_BACKGROUND6 = ndimage.generate_binary_structure(3, 1)


@dataclass(frozen=True)
class ThinACEImageResult:
    """RAS-oriented NIfTI result of the source ThinACE operation."""

    skeleton: nib.Nifti1Image


def thin_ace(skeleton: np.ndarray) -> np.ndarray:
    """Thin a positive ACE skeleton with distance-priority simple-point peeling.

    This follows the source's iteration structure, 18-connected foreground /
    6-connected background tests, endpoint retention, and zero image border.
    Output foreground is encoded as 254, matching `ThinACE.solve`.
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
    distance = _signed_distance_fmm(initial_phi, max_distance=4.0)

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
        heap: list[tuple[float, int, int, int]] = []
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
                heapq.heappush(heap, (float(distance[point]), k, j, i))
                found_simple = True
            else:
                labels[point] = 2
                list2.append(point)

        while heap:
            _, k, j, i = heapq.heappop(heap)
            point = (i, j, k)
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
            for dx, dy, dz in _OFFSETS_6:
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
    neighborhood = object_mask[
        point[0] - 1 : point[0] + 2,
        point[1] - 1 : point[1] + 2,
        point[2] - 1 : point[2] + 2,
    ].copy()
    neighborhood[1, 1, 1] = False
    foreground = neighborhood & _FOREGROUND18
    background = ~neighborhood & _BACKGROUND6
    background[1, 1, 1] = False
    return (
        _component_count(foreground, _FOREGROUND18),
        _component_count(background, _BACKGROUND6),
    )


def _component_count(mask: np.ndarray, structure: np.ndarray) -> int:
    _, count = ndimage.label(mask, structure=structure)
    return int(count)


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
