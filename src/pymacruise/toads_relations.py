"""TOADS relation-neighborhood maps from the selected VandyCRUISE source."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


# Source nested i/j/l loops visit these six faces in this order.
_FACES = ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0))
_INTERIOR = (slice(1, -1), slice(1, -1), slice(1, -1))
_EMPTY = np.int8(-1)


@dataclass(frozen=True)
class ToadsRelationList:
    """Source-style relation map and relation-membership table."""

    relation_map: NDArray[np.int16]
    relation_membership: NDArray[np.bool_]
    diffusion_iterations: int = 0

    @property
    def class_adjacency(self) -> NDArray[np.bool_]:
        """Return the selected source's first ``class_count`` union rows."""
        return self.relation_membership[: self.relation_membership.shape[1]]

    @property
    def class_count(self) -> int:
        return int(self.relation_membership.shape[1])


@dataclass(frozen=True)
class ToadsRelationState:
    """Composed relation, boundary-distance, and closest-label state."""

    relation: ToadsRelationList
    boundary_distances: NDArray[np.int8]
    closest_labels: NDArray[np.int8]
    boundary_diffusion_passes: int
    closest_label_passes: int


def toads_refresh_relation_state(
    segmentation: ArrayLike,
    memberships: ArrayLike,
    boundary_distances: ArrayLike,
    limit: int,
) -> ToadsRelationState:
    """Refresh the selected source's three related state products in order.

    ``AlgorithmMultipleToads`` invokes ``computeRelations``, then
    ``computeBoundaryDistances``, then ``computeClosestLabels`` after a
    segmentation update. This helper preserves that order, passes through the
    caller's carried boundary-distance state, and evaluates all products from
    the same validated segmentation snapshot. The boundary pass count follows
    the selected source's reset ``change`` flag (one pass when ``limit > 0``);
    closest labels execute the fixed ``5*limit - 1`` passes when positive.
    """
    labels = np.asarray(segmentation)
    membership_values = np.asarray(memberships)
    distances = np.asarray(boundary_distances)
    if labels.ndim != 3:
        raise ValueError("segmentation must be a 3-D array")
    if membership_values.ndim != 4 or membership_values.shape[:3] != labels.shape:
        raise ValueError("memberships must match the segmentation grid and have a class axis")
    if distances.ndim != 4 or distances.shape[:3] != labels.shape:
        raise ValueError("boundary_distances must match the segmentation grid and have a class axis")
    if distances.shape[3] != membership_values.shape[3]:
        raise ValueError("memberships and boundary_distances must have the same class count")
    if isinstance(limit, (bool, np.bool_)) or not isinstance(limit, (int, np.integer)) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")

    relation = toads_compute_relations(labels, membership_values, int(limit))
    updated_distances = toads_update_boundary_distances(labels, distances, int(limit))
    closest = toads_closest_labels(labels, int(limit))
    return ToadsRelationState(
        relation=relation,
        boundary_distances=updated_distances,
        closest_labels=closest,
        boundary_diffusion_passes=int(limit > 0),
        closest_label_passes=max(0, 5 * int(limit) - 1),
    )


def toads_distance_factor(boundary_distances: ArrayLike, marching_limit: float) -> NDArray[np.float32]:
    """Return the source ``1 / (1 + boundaryDist / limit)`` factor.

    The input is the signed-byte per-class field produced by
    ``computeBoundaryDistances``. Float32 division and addition follow Java's
    field types and operation order. The selected class accepts only a
    positive finite marching limit for this helper.
    """
    distances = np.asarray(boundary_distances)
    if not np.issubdtype(distances.dtype, np.integer) or np.issubdtype(distances.dtype, np.bool_):
        raise ValueError("boundary_distances must contain integer signed-byte values")
    bounds = np.iinfo(np.int8)
    if np.any(distances < bounds.min) or np.any(distances > bounds.max):
        raise ValueError("boundary_distances must fit Java signed-byte range")
    limit = float(marching_limit)
    if not np.isfinite(limit) or limit <= 0.0:
        raise ValueError("marching_limit must be finite and positive")
    distance_values = distances.astype(np.float32)
    limit_value = np.float32(limit)
    quotient = np.asarray(distance_values / limit_value, dtype=np.float32)
    denominator = np.asarray(np.float32(1.0) + quotient, dtype=np.float32)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return np.asarray(np.float32(1.0) / denominator, dtype=np.float32)


def toads_shape_factor(shape: ArrayLike, prior_factor: float) -> NDArray[np.float32]:
    """Return the selected ``min(1, 2*priorFactor*shape)`` prior factor."""
    values = np.asarray(shape, dtype=np.float32)
    prior = float(prior_factor)
    if not np.isfinite(values).all() or not np.isfinite(prior):
        raise ValueError("shape and prior_factor must be finite")
    scaled = np.asarray(np.float32(2.0) * np.float32(prior), dtype=np.float32)
    scaled = np.asarray(scaled * values, dtype=np.float32)
    # Match Numerics.min(a,b): it returns b when a < b is false.
    return np.asarray(np.where(np.float32(1.0) < scaled, np.float32(1.0), scaled), dtype=np.float32)


def toads_neighbor_factor(shape: ArrayLike, prior_factor: float) -> NDArray[np.float32]:
    """Return selected ``bounded(2*priorFactor*shape, 0.25, 1)`` semantics."""
    values = np.asarray(shape, dtype=np.float32)
    prior = float(prior_factor)
    if not np.isfinite(values).all() or not np.isfinite(prior):
        raise ValueError("shape and prior_factor must be finite")
    scaled = np.asarray(np.float32(2.0) * np.float32(prior), dtype=np.float32)
    scaled = np.asarray(scaled * values, dtype=np.float32)
    # Numerics.bounded is max(low, min(value, high)), using comparisons whose
    # NaN behavior differs from np.clip; finite inputs keep this source-exact.
    upper_bounded = np.where(scaled < np.float32(1.0), scaled, np.float32(1.0))
    lower_bounded = np.where(np.float32(0.25) > upper_bounded, np.float32(0.25), upper_bounded)
    return np.asarray(lower_bounded, dtype=np.float32)


def toads_membership_relation_tiers(
    closest_labels: ArrayLike,
    transformed_shapes: ArrayLike,
) -> tuple[NDArray[np.int8], NDArray[np.bool_]]:
    """Classify per-voxel classes as direct, neighbor, or prior support.

    The returned byte codes match `computeExactMemberships`: direct=2 for
    closest-label slots 0/1, neighbor=1 for slots 2/3, and prior=3 otherwise.
    The paired eligibility mask marks prior-tier classes as excluded only
    when the transformed shape value is exactly zero, matching the source's
    early zero assignment before computing that class's energy.
    """
    labels = np.asarray(closest_labels)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    if labels.ndim != 4 or labels.shape[-1] != 4 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("closest_labels must be an integer grid with four slots")
    if shapes.ndim != 4 or shapes.shape[:3] != labels.shape[:3] or shapes.shape[-1] < 1:
        raise ValueError("transformed_shapes must match the closest-label grid and have classes")
    if not np.isfinite(shapes).all():
        raise ValueError("transformed_shapes must be finite")
    if np.any(labels < -1) or np.any(labels >= shapes.shape[-1]):
        raise ValueError("closest_labels must contain only the -1 sentinel or valid class indices")

    direct = np.any(labels[..., :2, None] == np.arange(shapes.shape[-1]), axis=-2)
    neighbor = np.any(labels[..., 2:, None] == np.arange(shapes.shape[-1]), axis=-2) & ~direct
    relation_types = np.full(shapes.shape, 3, dtype=np.int8)
    relation_types[neighbor] = np.int8(1)
    relation_types[direct] = np.int8(2)
    eligible = ~((relation_types == 3) & (shapes == np.float32(0.0)))
    return relation_types, eligible


def toads_closest_labels(segmentation: ArrayLike, limit: int) -> NDArray[np.int8]:
    """Build the selected source's four-slot closest-label map.

    This implements ``MultipleToadSegmentation.computeClosestLabels``. Slot
    zero is the voxel's own class, slots one through three are distinct face
    neighbor labels, and missing entries are Java's ``EMPTY`` byte (-1).
    Missing slots diffuse synchronously for source passes ``1 <= t < 5*limit``.
    The selected method scans only interior voxels while seeding/diffusing;
    initialized border slots therefore remain unchanged. Operations are
    vectorized over voxels while preserving the Java face traversal order.
    """
    labels = np.asarray(segmentation)
    if labels.ndim != 3 or any(size == 0 for size in labels.shape):
        raise ValueError("segmentation must be a nonempty 3-D integer array")
    if not np.issubdtype(labels.dtype, np.integer) or np.issubdtype(labels.dtype, np.bool_):
        raise ValueError("segmentation labels must be integers")
    if any(size < 3 for size in labels.shape):
        raise ValueError("segmentation dimensions must be at least three voxels")
    if np.any(labels < 0) or np.any(labels > 127):
        raise ValueError("segmentation labels must fit the nonnegative Java byte class range")
    if isinstance(limit, (bool, np.bool_)) or not isinstance(limit, (int, np.integer)) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")

    class_labels = labels.astype(np.int8, copy=False)
    closest = np.full((*labels.shape, 4), _EMPTY, dtype=np.int8)
    closest[..., 0] = class_labels
    interior = closest[_INTERIOR]

    # Seed each interior voxel's distinct face labels with the selected loop order.
    for dx, dy, dz in _FACES:
        neighbor = class_labels[
            slice(1 + dx, labels.shape[0] - 1 + dx),
            slice(1 + dy, labels.shape[1] - 1 + dy),
            slice(1 + dz, labels.shape[2] - 1 + dz),
        ]
        slot1, slot2, slot3 = (interior[..., slot] for slot in (1, 2, 3))
        different = neighbor != interior[..., 0]
        fill1 = (slot1 == _EMPTY) & different
        fill2 = (~(slot1 == _EMPTY)) & (slot2 == _EMPTY) & (slot1 != neighbor) & different
        fill3 = (
            (~fill1)
            & (~fill2)
            & different
            & (slot3 == _EMPTY)
            & (slot1 != neighbor)
            & (slot2 != neighbor)
        )
        slot1[fill1] = neighbor[fill1]
        slot2[fill2] = neighbor[fill2]
        slot3[fill3] = neighbor[fill3]

    for _ in range(1, 5 * int(limit)):
        updated = closest.copy()
        target = updated[_INTERIOR]
        current = closest
        needs_slot1 = target[..., 1] == _EMPTY
        needs_slot2 = (target[..., 1] != _EMPTY) & (target[..., 2] == _EMPTY)
        needs_slot3 = (
            (target[..., 1] != _EMPTY)
            & (target[..., 2] != _EMPTY)
            & (target[..., 3] == _EMPTY)
        )
        for dx, dy, dz in _FACES:
            neighbor = current[
                slice(1 + dx, labels.shape[0] - 1 + dx),
                slice(1 + dy, labels.shape[1] - 1 + dy),
                slice(1 + dz, labels.shape[2] - 1 + dz),
                :,
            ]
            slot1 = target[..., 1]
            slot2 = target[..., 2]
            slot3 = target[..., 3]
            n1, n2, n3 = neighbor[..., 1], neighbor[..., 2], neighbor[..., 3]
            accept1 = needs_slot1 & (n1 != _EMPTY) & (target[..., 0] != n1)
            accept21 = (
                needs_slot2
                & (n1 != _EMPTY)
                & (target[..., 0] != n1)
                & (slot1 != n1)
            )
            accept22 = (
                needs_slot2
                & ~((n1 != _EMPTY) & (target[..., 0] != n1) & (slot1 != n1))
                & (n2 != _EMPTY)
                & (target[..., 0] != n2)
                & (slot1 != n2)
            )
            accept31 = (
                needs_slot3
                & (n1 != _EMPTY)
                & (target[..., 0] != n1)
                & (slot1 != n1)
                & (slot2 != n1)
            )
            accept32 = (
                needs_slot3
                & ~((n1 != _EMPTY) & (target[..., 0] != n1) & (slot1 != n1) & (slot2 != n1))
                & (n2 != _EMPTY)
                & (target[..., 0] != n2)
                & (slot1 != n2)
                & (slot2 != n2)
            )
            accept33 = (
                needs_slot3
                & ~(
                    ((n1 != _EMPTY) & (target[..., 0] != n1) & (slot1 != n1) & (slot2 != n1))
                    | ((n2 != _EMPTY) & (target[..., 0] != n2) & (slot1 != n2) & (slot2 != n2))
                )
                & (n3 != _EMPTY)
                & (target[..., 0] != n3)
                & (slot1 != n3)
                & (slot2 != n3)
            )
            slot1[accept1] = n1[accept1]
            slot2[accept21] = n1[accept21]
            slot2[accept22] = n2[accept22]
            slot3[accept31] = n1[accept31]
            slot3[accept32] = n2[accept32]
            slot3[accept33] = n3[accept33]
        closest[_INTERIOR] = target

    return closest


def toads_build_relation_list(
    segmentation: ArrayLike,
    class_count: int,
) -> ToadsRelationList:
    """Build source-ordered relation IDs for local six-face label sets.

    This ports ``MultipleToadSegmentation.buildRelationList`` from a fresh
    constructor-equivalent relation map. Singleton relations occupy IDs
    ``0..class_count-1``. Each distinct boundary set, containing the voxel's
    own class and its unique six-face neighbors, receives the next ID in the
    source x/y/z first-occurrence order; set equality is order-independent.
    The returned table contains those sets followed by the source's unioned
    adjacency rows for each singleton class. Borders retain their initialized
    class labels, as the Java method writes only interior voxels. A uint64
    bitset is used internally to vectorize boundary set collection.
    """
    labels = np.asarray(segmentation)
    if labels.ndim != 3 or any(size == 0 for size in labels.shape):
        raise ValueError("segmentation must be a nonempty 3-D integer array")
    if not np.issubdtype(labels.dtype, np.integer) or np.issubdtype(labels.dtype, np.bool_):
        raise ValueError("segmentation labels must be integers")
    if isinstance(class_count, (bool, np.bool_)) or not isinstance(class_count, (int, np.integer)):
        raise ValueError("class_count must be a positive integer")
    classes = int(class_count)
    if not 1 <= classes <= 63:
        raise ValueError("class_count must be between 1 and 63 for source class bitsets")
    if any(size < 3 for size in labels.shape):
        raise ValueError("segmentation dimensions must be at least three voxels")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation contains an out-of-range class index")

    class_labels = labels.astype(np.int64, copy=False)
    center = class_labels[_INTERIOR]
    relation_keys = np.left_shift(np.uint64(1), center.astype(np.uint64))
    is_boundary = np.zeros(center.shape, dtype=bool)
    for dx, dy, dz in _FACES:
        neighbor = class_labels[
            slice(1 + dx, labels.shape[0] - 1 + dx),
            slice(1 + dy, labels.shape[1] - 1 + dy),
            slice(1 + dz, labels.shape[2] - 1 + dz),
        ]
        differs = neighbor != center
        is_boundary |= differs
        bits = np.left_shift(np.uint64(1), neighbor.astype(np.uint64))
        relation_keys |= np.where(differs, bits, np.uint64(0))

    boundary_keys = relation_keys[is_boundary]
    if boundary_keys.size:
        unique_keys, first_indices = np.unique(boundary_keys, return_index=True)
        ordered_keys = unique_keys[np.argsort(first_indices)]
    else:
        ordered_keys = np.empty(0, dtype=np.uint64)
    if classes + ordered_keys.size > 500:
        raise ValueError("source relation table exceeds the 500-entry capacity")

    relation_map = class_labels.astype(np.int16, copy=True)
    relation_map_interior = relation_map[_INTERIOR]
    for relation_id, key in enumerate(ordered_keys, start=classes):
        relation_map_interior[is_boundary & (relation_keys == key)] = relation_id

    relation_membership = np.zeros((classes + ordered_keys.size, classes), dtype=bool)
    relation_membership[np.arange(classes), np.arange(classes)] = True
    for index, key in enumerate(ordered_keys, start=classes):
        members = np.array(
            [bool((int(key) >> class_index) & 1) for class_index in range(classes)],
            dtype=bool,
        )
        relation_membership[index] = members
        for class_index in np.flatnonzero(members):
            relation_membership[class_index] |= members

    return ToadsRelationList(relation_map, relation_membership)


def toads_compute_relations(
    segmentation: ArrayLike,
    memberships: ArrayLike,
    limit: int,
) -> ToadsRelationList:
    """Build and diffuse the selected exact-membership TOADS relation map.

    This composes the source `computeRelations` loop with the validated local
    relation table. Class-only interior entries may adopt a neighboring
    relation ID with the strictly greatest summed membership support. Face
    traversal and class accumulation match the Java ordering; each pass reads
    the prior map and writes synchronously. The source sets ``change`` while
    any simple class map remains and has an exclusive ``t < 10*limit`` cap.
    Approximate-membership mode is not covered and is intentionally rejected.
    """
    labels = np.asarray(segmentation)
    membership_values = np.asarray(memberships, dtype=np.float32)
    if membership_values.ndim != 4 or membership_values.shape[:3] != labels.shape:
        raise ValueError("memberships must match the segmentation grid and have a class axis")
    classes = membership_values.shape[-1]
    if classes < 1 or not np.isfinite(membership_values).all():
        raise ValueError("memberships must be finite and have at least one class")
    if isinstance(limit, (bool, np.bool_)) or not isinstance(limit, (int, np.integer)) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")

    relation = toads_build_relation_list(labels, classes)
    relation_map = relation.relation_map.copy()
    table = relation.relation_membership
    nx, ny, nz = labels.shape
    iterations = 0
    change = True
    for _ in range(1, 10 * int(limit)):
        if not change:
            break
        current = relation_map
        updated = relation_map.copy()
        target = updated[_INTERIOR]
        change = bool(np.any(current[_INTERIOR] < classes))
        active = target < classes
        best_support = np.zeros(active.shape, dtype=np.float32)
        best_relation = target.copy()
        for dx, dy, dz in _FACES:
            slices = (
                slice(1 + dx, nx - 1 + dx),
                slice(1 + dy, ny - 1 + dy),
                slice(1 + dz, nz - 1 + dz),
            )
            neighbor_relation = current[slices]
            candidate = neighbor_relation >= classes
            safe_relation = np.clip(neighbor_relation, 0, table.shape[0] - 1)
            related_classes = table[safe_relation]
            neighbor_memberships = membership_values[slices]
            support = np.zeros(active.shape, dtype=np.float32)
            for class_index in range(classes):
                term = np.where(
                    related_classes[..., class_index],
                    neighbor_memberships[..., class_index],
                    np.float32(0.0),
                ).astype(np.float32, copy=False)
                support = np.asarray(support + term, dtype=np.float32)
            choose = active & candidate & (support > best_support)
            best_support[choose] = support[choose]
            best_relation[choose] = neighbor_relation[choose]
        target[...] = best_relation
        relation_map[_INTERIOR] = target
        iterations += 1

    return ToadsRelationList(
        relation_map=relation_map,
        relation_membership=table,
        diffusion_iterations=iterations,
    )


def toads_update_boundary_distances(
    segmentation: ArrayLike,
    boundary_distances: ArrayLike,
    limit: int,
) -> NDArray[np.int8]:
    """Apply the selected source's class-wise TOADS boundary-distance update.

    The method seeds distance 1 on both class sides of every interior six-face
    transition and retains existing values. Its later diffusion excerpt is
    intentionally preserved: `change` is reset and never set true, while
    `Numerics.min(current_zero, neighbor_positive + 1)` remains zero. Thus the
    source performs one no-op diffusion pass when ``limit > 0`` and produces
    only interface seeds. This behavior is source evidence, not an improved
    distance transform.
    """
    labels = np.asarray(segmentation)
    distances = np.asarray(boundary_distances)
    if labels.ndim != 3 or any(size == 0 for size in labels.shape):
        raise ValueError("segmentation must be a nonempty 3-D integer array")
    if not np.issubdtype(labels.dtype, np.integer) or np.issubdtype(labels.dtype, np.bool_):
        raise ValueError("segmentation labels must be integers")
    if distances.ndim != 4 or distances.shape[:3] != labels.shape or distances.shape[3] < 1:
        raise ValueError("boundary_distances must match the segmentation grid and have a class axis")
    if not np.issubdtype(distances.dtype, np.integer) or np.issubdtype(distances.dtype, np.bool_):
        raise ValueError("boundary_distances must contain integer byte values")
    if any(size < 3 for size in labels.shape):
        raise ValueError("segmentation dimensions must be at least three voxels")
    if np.any(labels < 0) or np.any(labels >= distances.shape[3]):
        raise ValueError("segmentation contains an out-of-range class index")
    limits = np.iinfo(np.int8)
    if np.any(distances < limits.min) or np.any(distances > limits.max):
        raise ValueError("boundary_distances must fit Java signed-byte range")
    if isinstance(limit, (bool, np.bool_)) or not isinstance(limit, (int, np.integer)) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")

    class_labels = labels.astype(np.int64, copy=False)
    output = distances.astype(np.int8, copy=True)
    center_labels = class_labels[_INTERIOR]
    center_distances = output[_INTERIOR]
    for dx, dy, dz in _FACES:
        neighbor_slices = (
            slice(1 + dx, labels.shape[0] - 1 + dx),
            slice(1 + dy, labels.shape[1] - 1 + dy),
            slice(1 + dz, labels.shape[2] - 1 + dz),
        )
        neighbor_labels = class_labels[neighbor_slices]
        neighbor_distances = output[neighbor_slices]
        differs = center_labels != neighbor_labels
        for class_index in range(output.shape[3]):
            center_side = differs & (neighbor_labels == class_index)
            neighbor_side = differs & (center_labels == class_index)
            center_distances[..., class_index][center_side] = np.int8(1)
            neighbor_distances[..., class_index][neighbor_side] = np.int8(1)

    change = True
    for _ in range(1, 5 * int(limit)):
        if not change:
            break
        for class_index in range(output.shape[3]):
            change = False
            current = output[..., class_index]
            temp = current[_INTERIOR].copy()
            active = (center_labels != class_index) & (temp == 0)
            for dx, dy, dz in _FACES:
                neighbor_slices = (
                    slice(1 + dx, labels.shape[0] - 1 + dx),
                    slice(1 + dy, labels.shape[1] - 1 + dy),
                    slice(1 + dz, labels.shape[2] - 1 + dz),
                )
                neighbor_distance = current[neighbor_slices]
                positive = active & (neighbor_distance > 0)
                candidate = neighbor_distance.astype(np.int16) + 1
                temp[positive] = np.minimum(temp[positive].astype(np.int16), candidate[positive]).astype(np.int8)
            current[_INTERIOR] = temp

    return output
