"""LesionTOADS growing-frontier initialization.

This ports the selected VandyCRUISE ``LesionToadSegmentation`` methods
``initGrowingLabels``, ``bestMembership``, and ``maxMembershipRatio`` through
queue seeding only. It reuses the already source-validated ``BinaryHeap4D``
port in :mod:`pymacruise.toads_heap`; it does not implement queue events,
topology checks, or lesion segmentation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.toads_heap import ToadsPriorityHeap4D


_FACE_NEIGHBORS = ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0))
_ZERO = np.float32(1.0e-10)
_LESION_OBJECT_TYPES = {"lesion", "blackhole"}


@dataclass(frozen=True)
class LesionToadsGrowingInitialization:
    """Availability and min-heap entries after ``initGrowingLabels``."""

    available: NDArray[np.bool_]
    queue: ToadsPriorityHeap4D


def _lesion_score(
    ordinary_membership: np.float32,
    lesion_membership: np.float32,
    black_hole_membership: np.float32,
    *,
    use_lesion_max: bool,
) -> np.float32:
    if use_lesion_max:
        # VandyCRUISE Numerics.max(a,b,c) is max(a, max(b,c)), with strict
        # greater-than comparisons; values here are finite and this preserves
        # its float32 result and tie behavior.
        inner = lesion_membership if lesion_membership > black_hole_membership else black_hole_membership
        return ordinary_membership if ordinary_membership > inner else inner
    return np.float32(
        np.float32(ordinary_membership + lesion_membership) + black_hole_membership
    )


def initialize_lesion_toads_growing_frontier(
    skeleton: ArrayLike,
    memberships: ArrayLike,
    ordering: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    *,
    initial_available: ArrayLike | None = None,
    use_lesion_max: bool = False,
    empty_label: int = -1,
    zero: float = 1.0e-10,
) -> LesionToadsGrowingInitialization:
    """Seed the selected LesionTOADS growing queue from a class skeleton.

    Inputs use ``(x, y, z)`` voxels and ``(x, y, z, class)`` memberships.
    Interior voxels are reset available for every class iff their skeleton
    label is empty. Each available voxel that touches class ``k`` on a
    six-face is removed from availability for ``k`` and inserted with priority
    ``ordering + maxMembershipRatio``. Border availability is preserved from
    the input but must be false for the source's unguarded neighbor reads.

    The selected Java ``bestMembership`` has a source quirk: it compares the
    adjusted lesion score against its running maximum, but when it wins it
    stores the ordinary membership (not that adjusted score) as the next
    maximum. This function retains that operation order. Only source object
    types ``lesion`` and ``blackhole`` receive the special score; optimized
    lesion/black-hole types do not in these two methods.
    """
    skel = np.asarray(skeleton)
    mems = np.asarray(memberships, dtype=np.float32)
    orderings = np.asarray(ordering, dtype=np.float32)
    lesion_mems = np.asarray(lesion_membership, dtype=np.float32)
    black_hole_mems = np.asarray(black_hole_membership, dtype=np.float32)
    types = tuple(object_types)

    if skel.ndim != 3 or any(size < 3 for size in skel.shape):
        raise ValueError("skeleton must be a 3-D grid with every dimension at least 3")
    if not np.issubdtype(skel.dtype, np.integer):
        raise ValueError("skeleton labels must be integers")
    if mems.ndim != 4 or mems.shape[:3] != skel.shape or mems.shape[-1] < 2:
        raise ValueError("memberships must match skeleton and contain at least two classes")
    classes = int(mems.shape[-1])
    if classes > 127:
        raise ValueError("class indices must fit the Java signed-byte loop")
    if len(types) != classes or not all(isinstance(value, str) and value for value in types):
        raise ValueError("one nonempty object type is required for every class")
    if orderings.shape != skel.shape:
        raise ValueError("ordering must match the skeleton grid")
    for name, values in (
        ("memberships", mems),
        ("ordering", orderings),
        ("lesion_membership", lesion_mems),
        ("black_hole_membership", black_hole_mems),
    ):
        if not np.isfinite(values).all():
            raise ValueError(f"{name} must contain only finite values")
    if lesion_mems.shape != skel.shape or black_hole_mems.shape != skel.shape:
        raise ValueError("lesion and black-hole memberships must match skeleton")
    if not isinstance(empty_label, (int, np.integer)):
        raise ValueError("empty_label must be an integer")
    if np.any((skel != empty_label) & ((skel < 0) | (skel >= classes))):
        raise ValueError("skeleton labels must be valid class indices or empty_label")
    zero_value = np.float32(zero)
    if not np.isfinite(zero_value) or zero_value < np.float32(0.0):
        raise ValueError("zero must be finite and nonnegative")

    if initial_available is None:
        available = np.zeros(skel.shape + (classes,), dtype=np.bool_)
    else:
        supplied = np.asarray(initial_available)
        if supplied.shape != skel.shape + (classes,) or supplied.dtype != np.bool_:
            raise ValueError("initial_available must be boolean and match skeleton and classes")
        available = supplied.copy()
    border = np.zeros(skel.shape, dtype=np.bool_)
    for axis in range(3):
        edge = [slice(None)] * 3
        edge[axis] = 0
        border[tuple(edge)] = True
        edge[axis] = -1
        border[tuple(edge)] = True
    if np.any(available[border]):
        raise ValueError("initial border availability must be false for source-safe neighbor reads")

    nx, ny, nz = skel.shape
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                available[x, y, z, :] = int(skel[x, y, z]) == empty_label

    queue = ToadsPriorityHeap4D(nx * ny + ny * nz + nz * nx, "min")
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for label in range(classes):
                    if not available[x, y, z, label]:
                        continue
                    if not any(
                        int(skel[x + dx, y + dy, z + dz]) == label
                        for dx, dy, dz in _FACE_NEIGHBORS
                    ):
                        continue
                    available[x, y, z, label] = False

                    best = 0
                    best_running_value = np.float32(0.0)
                    for candidate in range(classes):
                        if types[candidate] in _LESION_OBJECT_TYPES:
                            candidate_score = _lesion_score(
                                mems[x, y, z, candidate],
                                lesion_mems[x, y, z],
                                black_hole_mems[x, y, z],
                                use_lesion_max=use_lesion_max,
                            )
                        else:
                            candidate_score = mems[x, y, z, candidate]
                        if candidate_score > best_running_value:
                            best = candidate
                            # Deliberately mirrors source: do not assign
                            # candidate_score here; Java stores mems[...,m].
                            best_running_value = mems[x, y, z, candidate]

                    if label == best:
                        ratio = np.float32(1.0)
                    else:
                        if types[label] in _LESION_OBJECT_TYPES:
                            current = _lesion_score(
                                mems[x, y, z, label],
                                lesion_mems[x, y, z],
                                black_hole_mems[x, y, z],
                                use_lesion_max=use_lesion_max,
                            )
                        else:
                            current = mems[x, y, z, label]
                        if types[best] in _LESION_OBJECT_TYPES:
                            maximum = _lesion_score(
                                mems[x, y, z, best],
                                lesion_mems[x, y, z],
                                black_hole_mems[x, y, z],
                                use_lesion_max=use_lesion_max,
                            )
                        else:
                            maximum = mems[x, y, z, best]
                        numerator = np.float32(maximum + zero_value)
                        denominator = np.float32(current + zero_value)
                        ratio = np.float32(numerator / denominator)
                    priority = np.float32(orderings[x, y, z] + ratio)
                    queue.add(priority, x, y, z, label)

    return LesionToadsGrowingInitialization(available=available, queue=queue)
