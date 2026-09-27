"""Native Python port of VandyCRUISE's 3D (6,18) topology correction.

Source anchors are the topology implementation by Pierre-Louis Bazin and
Blake Lucas in VandyCRUISE's `TopologyCorrection`, `TopologyPropagation`, and
`NeighborConnectivity`, plus the Vanderbilt MaCRUISE adapters. Their supplied
notices and attribution inventory are preserved under `licenses/upstream/` and
`license_info.md`.
"""

from __future__ import annotations

from dataclasses import dataclass
import heapq
from itertools import count

import nibabel as nib
import numpy as np
from scipy import ndimage

from pymacruise.image_io import normalize_image_to_ras

_FRONTIER = 0
_MASK = 10
_CRITICAL = 20
_BACKGROUND = 30
_BOUNDARY = 40
_TEST = 70
_OBJECT = 80

_FOREGROUND_6 = ndimage.generate_binary_structure(3, 1)
_BACKGROUND_18 = ndimage.generate_binary_structure(3, 2)


@dataclass(frozen=True)
class TopologyCorrectionResult:
    """Topology-corrected scalar field and optional atlas-refined product."""

    corrected: nib.Nifti1Image
    atlas_refined: nib.Nifti1Image | None


def topology_numbers_618(neighborhood: np.ndarray) -> tuple[int, int]:
    """Return local foreground-6 and background-18 component counts.

    ``neighborhood`` is a boolean 3x3x3 block around a candidate voxel. Its
    center is excluded from both counts, as in VandyCRUISE's
    `NeighborConnectivity` implementation.
    """
    values = np.asarray(neighborhood, dtype=bool)
    if values.shape != (3, 3, 3):
        raise ValueError("A topology neighborhood must have shape (3, 3, 3)")
    foreground = values.copy()
    background = ~values
    foreground[1, 1, 1] = False
    background[1, 1, 1] = False
    foreground[0, 0, 0] = False
    foreground[0, 0, 2] = False
    foreground[0, 2, 0] = False
    foreground[0, 2, 2] = False
    foreground[2, 0, 0] = False
    foreground[2, 0, 2] = False
    foreground[2, 2, 0] = False
    foreground[2, 2, 2] = False
    return (
        _count_components(foreground, _FOREGROUND_6),
        _count_components(background, _BACKGROUND_18),
    )


def topology_correct(
    volume: nib.spatialimages.SpatialImage,
    paint_mask: nib.spatialimages.SpatialImage,
    *,
    paint_threshold: float = 0.5,
    lowest_level: float = 0.0,
    highest_level: float = 1.0,
    regularize: float = 1e-4,
    multi_atlas_wm: nib.spatialimages.SpatialImage | None = None,
) -> TopologyCorrectionResult:
    """Apply the VandyCRUISE fast-marching (6,18) scalar correction.

    The supplied paint mask initializes the object; every positive paint voxel
    is set to the normalized highest level. Propagation proceeds downward in
    intensity order and admits a candidate only when adding it preserves one
    foreground 6-component and one background 18-component in its 3x3x3
    neighborhood. This follows the selected MaCRUISE JIST adapter defaults.

    If ``multi_atlas_wm`` is supplied, also reproduce V3's subsequent
    `Refine_Topo_by_MA` operation: first apply the source's one-pass 3D
    six-neighbor integer minimum filter (with zero-valued image borders), then
    set corrected values below 200 to 255 wherever that eroded WM is positive.
    """
    source = _as_ras_3d(volume, "topology volume")
    paint = _as_ras_3d(paint_mask, "topology paint mask")
    _require_same_grid(source, paint, "topology volume", "paint mask")
    if multi_atlas_wm is not None:
        atlas = _as_ras_3d(multi_atlas_wm, "multi-atlas WM mask")
        _require_same_grid(source, atlas, "topology volume", "multi-atlas WM mask")
        eroded_atlas = _erode_multi_atlas_wm(np.asarray(atlas.dataobj))
    else:
        eroded_atlas = None

    if not 0.0 <= paint_threshold or not np.isfinite(paint_threshold):
        raise ValueError("paint_threshold must be finite and non-negative")
    if not np.isfinite(lowest_level) or not np.isfinite(highest_level):
        raise ValueError("Topology levels must be finite")
    if highest_level <= lowest_level:
        raise ValueError("highest_level must be greater than lowest_level")
    if not np.isfinite(regularize) or regularize < 0:
        raise ValueError("regularize must be finite and non-negative")

    image = np.asarray(source.dataobj, dtype=np.float32)
    mask = np.asarray(paint.dataobj) > paint_threshold
    if not np.isfinite(image).all():
        raise ValueError("Topology volume contains NaN or infinite voxel values")
    if not mask.any():
        raise ValueError("Topology paint mask contains no voxels above threshold")

    # AlgorithmTopologyCorrection adds a one-voxel zero frame before
    # computing its intensity range.
    minimum = min(float(image.min()), 0.0)
    maximum = max(float(image.max()), 0.0)
    if maximum <= minimum:
        # The legacy min/max-scaled regularization is zero for a constant image.
        corrected_data = image.copy()
    else:
        scaled_low = minimum + lowest_level * (maximum - minimum)
        scaled_high = maximum - (1.0 - highest_level) * (maximum - minimum)
        image = image.copy()
        image[mask] = scaled_high
        image[(image > scaled_high) | (image < scaled_low)] = scaled_low
        corrected_data = _downward_topology_propagation(
            image,
            mask,
            min_distance=regularize * (maximum - minimum),
            lowest_level=scaled_low,
        )

    corrected = _make_output(corrected_data, source)
    atlas_refined = None
    if eroded_atlas is not None:
        refined_data = corrected_data.copy()
        refined_data[(refined_data < 200.0) & (eroded_atlas > 0)] = 255.0
        atlas_refined = _make_output(refined_data, source)
    return TopologyCorrectionResult(corrected, atlas_refined)


def _downward_topology_propagation(
    image: np.ndarray,
    paint: np.ndarray,
    *,
    min_distance: float,
    lowest_level: float,
) -> np.ndarray:
    """Port `initDownwardPaintLabels` + exact downward smoothing (3D)."""
    nx, ny, nz = image.shape
    padded = np.zeros((nx + 2, ny + 2, nz + 2), dtype=np.float32)
    padded[1:-1, 1:-1, 1:-1] = image
    paint_padded = np.zeros(padded.shape, dtype=bool)
    paint_padded[1:-1, 1:-1, 1:-1] = paint
    object_mask = np.zeros(padded.shape, dtype=bool)
    object_mask[1:-1, 1:-1, 1:-1] = True

    labels = np.full(padded.shape, _BACKGROUND, dtype=np.uint8)
    distance = np.zeros(padded.shape, dtype=np.float32)
    labels[[0, -1], :, :] = _FRONTIER
    labels[:, [0, -1], :] = _FRONTIER
    labels[:, :, [0, -1]] = _FRONTIER
    labels[~object_mask] = _FRONTIER
    paint_padded &= object_mask

    heap: list[tuple[float, int, tuple[int, int, int]]] = []
    tie = count()
    labels[~object_mask] = _FRONTIER
    labels[object_mask & ~paint_padded] = _BACKGROUND
    labels[paint_padded] = _OBJECT
    distance[paint_padded] = padded[paint_padded]

    # Initial boundary is the 6-neighborhood of the paint mask.
    candidates = ndimage.binary_dilation(paint_padded, structure=_FOREGROUND_6)
    candidates &= object_mask & ~paint_padded
    labels[candidates] = _BOUNDARY
    distance[candidates] = padded[candidates]
    for xyz in np.argwhere(candidates):
        point = tuple(int(v) for v in xyz)
        heapq.heappush(heap, (-float(distance[point]), next(tie), point))

    last_value = 0.0
    while heap:
        neg_value, _, point = heapq.heappop(heap)
        value = -neg_value
        last_value = value
        x, y, z = point
        lower = _lowest_object_neighbor(distance, labels, point, value)
        if lower < value:
            heapq.heappush(heap, (-lower, next(tie), point))
            continue

        # The source adapter enables threshold stopping: downward propagation
        # does not promote points whose value has reached the global minimum.
        if value <= lowest_level:
            continue

        labels[point] = _TEST
        if _is_regular_object_addition(labels, point):
            labels[point] = _OBJECT
            distance[point] = value
            last_value = value
            for dx, dy, dz in _NEIGHBORS_18:
                neighbor = (x + dx, y + dy, z + dz)
                if labels[neighbor] == _CRITICAL:
                    labels[neighbor] = _BOUNDARY
                    heapq.heappush(heap, (-value, next(tie), neighbor))
                if (dx, dy, dz) in _FACE_OFFSETS and labels[neighbor] == _BACKGROUND:
                    labels[neighbor] = _BOUNDARY
                    proposed = min(float(padded[neighbor]), value - min_distance)
                    heapq.heappush(heap, (-proposed, next(tie), neighbor))
        else:
            labels[point] = _CRITICAL

    internal = (labels == _BACKGROUND) | (labels == _CRITICAL)
    distance[internal] = last_value
    return distance[1:-1, 1:-1, 1:-1].copy()


def _lowest_object_neighbor(
    distance: np.ndarray,
    labels: np.ndarray,
    point: tuple[int, int, int],
    value: float,
) -> float:
    x, y, z = point
    values = [
        float(distance[x + dx, y + dy, z + dz])
        for dx, dy, dz in _FACE_OFFSETS
        if labels[x + dx, y + dy, z + dz] >= _OBJECT
    ]
    return min(values, default=value)


def _is_regular_object_addition(
    labels: np.ndarray,
    point: tuple[int, int, int],
) -> bool:
    x, y, z = point
    neighborhood = labels[x - 1 : x + 2, y - 1 : y + 2, z - 1 : z + 2]
    candidate = neighborhood >= _TEST
    foreground, background = topology_numbers_618(candidate)
    return foreground == 1 and background == 1


def _count_components(mask: np.ndarray, structure: np.ndarray) -> int:
    labels, count_components = ndimage.label(mask, structure=structure)
    del labels
    return int(count_components)


def _erode_multi_atlas_wm(data: np.ndarray) -> np.ndarray:
    """Match MASImorph's default 3D center-plus-six minimum filter."""
    integer_data = np.asarray(data).astype(np.int32)
    return ndimage.minimum_filter(
        integer_data,
        footprint=_FOREGROUND_6,
        mode="constant",
        cval=0,
    )


_FACE_OFFSETS = (
    (-1, 0, 0),
    (1, 0, 0),
    (0, -1, 0),
    (0, 1, 0),
    (0, 0, -1),
    (0, 0, 1),
)

_NEIGHBORS_18 = tuple(
    (dx, dy, dz)
    for dx in (-1, 0, 1)
    for dy in (-1, 0, 1)
    for dz in (-1, 0, 1)
    if 0 < abs(dx) + abs(dy) + abs(dz) <= 2
)


def _as_ras_3d(
    image: nib.spatialimages.SpatialImage,
    name: str,
) -> nib.spatialimages.SpatialImage:
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError(f"{name} must be a NIfTI image")
    result = normalize_image_to_ras(image, source_name=name)
    if len(result.shape) != 3:
        raise ValueError(f"{name} must be a 3D scalar image")
    return result


def _require_same_grid(
    image_a: nib.spatialimages.SpatialImage,
    image_b: nib.spatialimages.SpatialImage,
    name_a: str,
    name_b: str,
) -> None:
    if image_a.shape != image_b.shape or not np.allclose(
        image_a.affine, image_b.affine, rtol=0.0, atol=1e-5
    ):
        raise ValueError(f"{name_a} and {name_b} must have the same physical grid")


def _make_output(data: np.ndarray, reference: nib.spatialimages.SpatialImage) -> nib.Nifti1Image:
    header = reference.header.copy()
    header.set_data_dtype(data.dtype)
    header.set_slope_inter(None, None)
    return nib.Nifti1Image(data, reference.affine, header=header)
