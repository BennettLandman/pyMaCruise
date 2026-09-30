"""Selected LesionTOADS relation-list construction and diffusion."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


# buildRelationList visits offsets in nested i/j/l order, retaining only the
# six face offsets. computeRelations has its own explicit if/else face order.
_BUILD_FACES = ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0))
_DIFFUSION_FACES = ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1))
_INTERIOR = (slice(1, -1), slice(1, -1), slice(1, -1))
_MAX_RELATIONS = 500


@dataclass(frozen=True)
class LesionRelationList:
    """The selected Java short relation map and boolean relation table."""

    relation_map: NDArray[np.int16]
    relation_membership: NDArray[np.bool_]
    diffusion_iterations: int = 0

    @property
    def relation_count(self) -> int:
        return int(self.relation_membership.shape[0])

    @property
    def class_adjacency(self) -> NDArray[np.bool_]:
        """Return the class rows after the source's boundary-relation union."""
        return self.relation_membership[: self.relation_membership.shape[1]]


def lesion_build_relation_list(segmentation: ArrayLike, class_count: int) -> LesionRelationList:
    """Build LesionTOADS relation IDs from selected ``buildRelationList``.

    Singleton class relations occupy the first ``class_count`` IDs. Interior
    boundary relation sets are assigned IDs at their first x/y/z occurrence;
    set matching ignores member order, while the six-face scan determines
    insertion order for newly encountered sets. Border voxels retain the
    constructor-initialized class labels. Relation rows are followed by class
    adjacency rows unioned from each multi-class boundary relation.
    """
    labels, classes = _validate(segmentation, class_count)
    if classes > _MAX_RELATIONS:
        raise ValueError(f"class_count exceeds source relation capacity ({_MAX_RELATIONS})")

    relation_map = labels.astype(np.int16, copy=True)
    relation_sets: list[frozenset[int]] = [frozenset((label,)) for label in range(classes)]
    relation_ids = {members: index for index, members in enumerate(relation_sets)}

    nx, ny, nz = labels.shape
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                center = int(labels[x, y, z])
                boundary = set()
                for dx, dy, dz in _BUILD_FACES:
                    neighbor = int(labels[x + dx, y + dy, z + dz])
                    if neighbor != center:
                        boundary.add(neighbor)
                if not boundary:
                    relation_map[x, y, z] = center
                    continue
                members = frozenset((*boundary, center))
                relation_id = relation_ids.get(members)
                if relation_id is None:
                    if len(relation_sets) >= _MAX_RELATIONS:
                        raise ValueError(
                            f"relation count exceeds source relation capacity ({_MAX_RELATIONS})"
                        )
                    relation_id = len(relation_sets)
                    relation_ids[members] = relation_id
                    relation_sets.append(members)
                relation_map[x, y, z] = relation_id

    relation_membership = np.zeros((len(relation_sets), classes), dtype=np.bool_)
    for relation_id, members in enumerate(relation_sets):
        for label in members:
            relation_membership[relation_id, label] = True
    # The source expands each singleton class row with every class found in a
    # multi-class relation containing that class.
    for class_index in range(classes):
        related_rows = relation_membership[classes:, class_index]
        if np.any(related_rows):
            relation_membership[class_index] |= np.any(
                relation_membership[classes:][related_rows], axis=0
            )

    return LesionRelationList(relation_map, relation_membership)


def lesion_compute_relations(
    segmentation: ArrayLike, class_count: int, limit: int
) -> LesionRelationList:
    """Build and synchronously diffuse relation IDs per selected Java methods.

    ``computeRelations`` runs ``buildRelationList`` first, then performs
    source passes ``t = 1 .. 5*limit-1``. Each class-only interior voxel takes
    the first neighboring boundary relation in X-, X+, Y-, Y+, Z-, Z+ order.
    The complete pass reads the previous relation map and is copied back only
    after all interiors have been calculated. Relation membership rows do not
    change during diffusion. This is a selected-method port, not a full
    LesionTOADS execution.
    """
    labels, classes = _validate(segmentation, class_count)
    if isinstance(limit, (bool, np.bool_)) or not isinstance(limit, (int, np.integer)) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")
    relation = lesion_build_relation_list(labels, classes)
    relation_map = relation.relation_map.copy()
    iterations = 0
    nx, ny, nz = labels.shape
    for _ in range(1, 5 * int(limit)):
        updated = relation_map.copy()
        for x in range(1, nx - 1):
            for y in range(1, ny - 1):
                for z in range(1, nz - 1):
                    if relation_map[x, y, z] < classes:
                        for dx, dy, dz in _DIFFUSION_FACES:
                            neighbor = int(relation_map[x + dx, y + dy, z + dz])
                            if neighbor >= classes:
                                updated[x, y, z] = neighbor
                                break
        relation_map = updated
        iterations += 1
    return LesionRelationList(
        relation_map, relation.relation_membership, diffusion_iterations=iterations
    )


def _validate(segmentation: ArrayLike, class_count: int) -> tuple[NDArray, int]:
    labels = np.asarray(segmentation)
    if labels.ndim != 3 or any(size < 3 for size in labels.shape):
        raise ValueError("segmentation must be a 3-D array with every size at least three")
    if not np.issubdtype(labels.dtype, np.integer) or np.issubdtype(labels.dtype, np.bool_):
        raise ValueError("segmentation labels must be integers")
    if isinstance(class_count, (bool, np.bool_)) or not isinstance(class_count, (int, np.integer)):
        raise ValueError("class_count must be an integer")
    classes = int(class_count)
    if not 1 <= classes <= 128:
        raise ValueError("class_count must fit LesionTOADS nonnegative signed-byte class labels")
    if np.any(labels < 0) or np.any(labels >= classes):
        raise ValueError("segmentation contains an out-of-range class index")
    return labels, classes
