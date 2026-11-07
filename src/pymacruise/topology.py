"""Native Python port of VandyCRUISE's 3D (6,18) topology correction.

Method context: the scalar-field, critical-point and fast-marching approach is
described by Bazin & Pham (2007); see ``docs/papers/topology-fast-marching.md``
for a code-to-paper map. This is not Han et al.'s graph-based GTCA: no graph
construction, cut direction, or multiscale graph optimization is performed
here. The concrete connectivity and numeric defaults below come from the
selected VandyCRUISE adapter, rather than being universal paper parameters.

Source anchors are the topology implementation by Pierre-Louis Bazin and
Blake Lucas in VandyCRUISE's `TopologyCorrection`, `TopologyPropagation`, and
`NeighborConnectivity`, plus the Vanderbilt MaCRUISE adapters. Their supplied
notices and attribution inventory are preserved under `licenses/upstream/` and
`license_info.md`.

Validation scope: the first direct harness comparison to selected Java
`AlgorithmTopologyCorrection` now matches all 729 voxels exactly, including
source max-tree tie order, threshold-stop fill, and float32 arithmetic. P25
remains in progress until a broader synthetic matrix exercises boundary,
plateau/tie, mask, and critical-point cases. See
`docs/topology-correction.md` for comparison evidence and limits.
"""

from __future__ import annotations

from dataclasses import dataclass

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
# Adapter-selected complementary digital connectivities. The paper describes
# enforcing a requested topology; the (6,18) pair is the VandyCRUISE/MaCRUISE
# configuration reproduced here, not a fixed parameter of every application.


class _SourceMaxTree:
    """Max-priority tree matching VandyCRUISE ``BinaryTree`` tie ordering."""

    def __init__(self) -> None:
        self._values = [np.float32(0.0)]
        self._points: list[tuple[int, int, int] | None] = [None]

    def __bool__(self) -> bool:
        return len(self._values) > 1

    def add(self, value: float, point: tuple[int, int, int]) -> None:
        node = len(self._values)
        parent = node // 2
        source_value = np.float32(value)
        self._values.append(np.float32(0.0))
        self._points.append(None)
        while parent > 0 and source_value > self._values[parent]:
            self._values[node] = self._values[parent]
            self._points[node] = self._points[parent]
            node = parent
            parent = node // 2
        self._values[node] = source_value
        self._points[node] = point

    def pop(self) -> tuple[np.float32, tuple[int, int, int]]:
        if not self:
            raise IndexError("pop from empty source max tree")
        first_value = self._values[1]
        first_point = self._points[1]
        assert first_point is not None
        last_value = self._values.pop()
        last_point = self._points.pop()
        size = len(self._values) - 1
        if size == 0:
            return first_value, first_point

        node = 1
        child = 2
        # BinaryTree.removeFirst uses strict comparisons and only considers
        # children before the final node being sifted down. Equal children
        # choose the left branch.
        while child <= size:
            if child + 1 <= size:
                if last_value < self._values[child] or last_value < self._values[child + 1]:
                    if self._values[child] < self._values[child + 1]:
                        promoted = child + 1
                    else:
                        promoted = child
                    self._values[node] = self._values[promoted]
                    self._points[node] = self._points[promoted]
                    node = promoted
                    child = node * 2
                else:
                    break
            else:
                if last_value < self._values[child]:
                    self._values[node] = self._values[child]
                    self._points[node] = self._points[child]
                    node = child
                    break
                break
        self._values[node] = last_value
        self._points[node] = last_point
        return first_value, first_point


@dataclass(frozen=True)
class TopologyCorrectionResult:
    """Topology-corrected scalar field and optional atlas-refined product."""

    corrected: nib.Nifti1Image
    atlas_refined: nib.Nifti1Image | None


def topology_numbers_618(neighborhood: np.ndarray) -> tuple[int, int]:
    """Return VandyCRUISE's local foreground-6 and background-18 counts.

    ``neighborhood`` is a boolean 3x3x3 block around a candidate voxel. Its
    center is excluded from both counts, as in VandyCRUISE's
    `NeighborConnectivity` implementation. Foreground components count only
    when they touch a face-neighbor of the center; isolated edge components
    are excluded by the source's `connected618Neighborhood3D`. These local
    counts are the discrete regularity test used by propagation. They are not
    Han et al.'s graph-based GTCA criterion.
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
    foreground_labels, _ = ndimage.label(foreground, structure=_FOREGROUND_6)
    face_component_labels = {
        int(foreground_labels[point])
        for point in (
            (0, 1, 1),
            (2, 1, 1),
            (1, 0, 1),
            (1, 2, 1),
            (1, 1, 0),
            (1, 1, 2),
        )
        if foreground_labels[point] > 0
    }
    return (
        len(face_component_labels),
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
    propagation_direction: str = "object->background",
    source_axis_order: str = "RAS",
) -> TopologyCorrectionResult:
    """Apply the current Python (6,18) scalar correction implementation.

    The supplied paint mask initializes the object; every positive paint voxel
    is set to the normalized highest level. Propagation proceeds downward in
    intensity order and admits a candidate only when adding it preserves one
    foreground 6-component and one background 18-component in its 3x3x3
    neighborhood. This follows the selected MaCRUISE JIST adapter defaults.
    In the paper's terms, the input is a scalar/implicit field and propagation
    preserves the requested topology across its level sets. The defaults
    `paint_threshold`, `lowest_level`, `highest_level`, and `regularize`
    reproduce wrapper choices; they are configurable adapter settings, not
    constants prescribed by Bazin & Pham. The optional multi-atlas WM
    refinement below is a separate V3 wrapper post-processing operation, not
    part of that paper's fast-marching method.

    Five selected-source synthetic fixtures now match exactly, including
    every queue event and critical-candidate count. This validates the tested
    scalar paint-mask cases, not whole-workflow or clinical equivalence; see
    the P25 report for scope.

    If ``multi_atlas_wm`` is supplied, also reproduce V3's subsequent
    `Refine_Topo_by_MA` operation: first apply the source's one-pass 3D
    six-neighbor integer minimum filter (with zero-valued image borders), then
    set corrected values below 200 to 255 wherever that eroded WM is positive.

    ``background->object`` selects the v3.2 intensity-initialized upward
    branch. Like the JIST wrapper, it ignores the paint values and seeds the
    expanded bounding box. ``source_axis_order=LPS`` preserves source queue
    tie order internally while returning the original public RAS grid. The
    existing paint-mask branch remains the default for Dura/shared callers.
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

    if propagation_direction not in ("object->background", "background->object"):
        raise ValueError("unsupported topology propagation direction")
    if source_axis_order not in ("RAS", "LPS"):
        raise ValueError("topology source_axis_order must be RAS or LPS")
    image = np.asarray(source.dataobj, dtype=np.float32)
    mask = np.asarray(paint.dataobj) > paint_threshold
    if not np.isfinite(image).all():
        raise ValueError("Topology volume contains NaN or infinite voxel values")
    if not mask.any() and propagation_direction == "object->background":
        raise ValueError("Topology paint mask contains no voxels above threshold")

    # AlgorithmTopologyCorrection adds a one-voxel zero frame before
    # computing its intensity range.
    minimum = min(float(image.min()), 0.0)
    maximum = max(float(image.max()), 0.0)
    if maximum <= minimum:
        # The legacy min/max-scaled regularization is zero for a constant image.
        corrected_data = image.copy()
    else:
        source_minimum = np.float32(minimum)
        source_maximum = np.float32(maximum)
        source_range = np.float32(source_maximum - source_minimum)
        source_low_level = np.float32(lowest_level)
        source_high_level = np.float32(highest_level)
        scaled_low = np.float32(
            source_minimum + source_low_level * source_range
        )
        scaled_high = np.float32(
            source_maximum
            - (np.float32(1.0) - source_high_level) * source_range
        )
        image = image.copy()
        if propagation_direction == "object->background":
            image[mask] = scaled_high
        image[(image > scaled_high) | (image < scaled_low)] = scaled_low
        if propagation_direction == "background->object":
            from .topology_upward import upward_exact_smoothing

            if any(size >= 32766 for size in image.shape):
                raise ValueError("topology indexing dimensions must be below 32766")
            values = image[::-1, ::-1, :] if source_axis_order == "LPS" else image
            corrected_data = upward_exact_smoothing(
                np.ascontiguousarray(values), scaled_low, scaled_high,
                np.float32(np.float32(regularize) * source_range),
            )
            if source_axis_order == "LPS":
                corrected_data = corrected_data[::-1, ::-1, :].copy()
        else:
            if source_axis_order != "RAS":
                raise ValueError("LPS topology order requires background->object propagation")
            corrected_data = _downward_topology_propagation(
                image, mask,
                min_distance=float(np.float32(np.float32(regularize) * source_range)),
                lowest_level=float(scaled_low),
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
    """Port `initDownwardPaintLabels` + exact downward smoothing (3D).

    This queue-driven descending front is the implementation-level
    fast-marching/topology-propagation mapping: candidate scalar values
    determine processing order, and an attempted object addition is accepted
    only when the local foreground/background component test remains regular.
    Rejected candidates are marked critical and are not promoted into the
    object. The one-voxel zero frame and minimum-level stopping follow the
    source adapter's boundary and threshold behavior; neither should be read
    as a paper-wide universal convention.
    """
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

    heap = _SourceMaxTree()
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
        heap.add(float(distance[point]), point)

    while heap:
        source_value, point = heap.pop()
        value = float(source_value)
        x, y, z = point
        lower = _lowest_object_neighbor(distance, labels, point, value)
        if lower < value:
            heap.add(lower, point)
            continue

        # The source adapter enables threshold stopping: downward propagation
        # does not promote points whose value has reached the global minimum.
        if value <= lowest_level:
            continue

        labels[point] = _TEST
        if _is_regular_object_addition(labels, point):
            labels[point] = _OBJECT
            distance[point] = value
            for dx, dy, dz in _NEIGHBORS_18:
                neighbor = (x + dx, y + dy, z + dz)
                if labels[neighbor] == _CRITICAL:
                    labels[neighbor] = _BOUNDARY
                    heap.add(value, neighbor)
                if (dx, dy, dz) in _FACE_OFFSETS and labels[neighbor] == _BACKGROUND:
                    labels[neighbor] = _BOUNDARY
                    proposed = np.float32(
                        min(float(padded[neighbor]), float(np.float32(value - min_distance)))
                    )
                    heap.add(float(proposed), neighbor)
        else:
            labels[point] = _CRITICAL

    internal = (labels == _BACKGROUND) | (labels == _CRITICAL)
    # The adapter enables TopologyPropagation.thresholdStop. Its final pass
    # fills unprocessed background and critical voxels with upLevel (the
    # selected normalized lowest level), not the last queue priority.
    distance[internal] = np.float32(lowest_level)
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
    """Apply the selected 3D digital-topology regularity criterion.

    The two component counts instantiate the adapter's (6,18) connectivity
    pair for a local addition. This is a local topology test within the
    fast-marching method; it is not the graph/cut optimization in Han et al.'s
    GTCA.
    """
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
