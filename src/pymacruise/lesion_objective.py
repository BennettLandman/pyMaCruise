"""Selected LesionTOADS objective/energy reduction over updated maps."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LesionObjectiveReduction:
    """Selected scalar objective and its source-order running total map."""

    total: np.float32
    running_total: np.ndarray


def reduce_lesion_toads_objective(
    normalized_class_memberships: np.ndarray,
    class_energies: np.ndarray,
    denominator: np.ndarray,
    lesion_membership: np.ndarray,
    black_hole_membership: np.ndarray,
    lesion_energy: np.ndarray,
    black_hole_energy: np.ndarray,
) -> LesionObjectiveReduction:
    """Accumulate the active selected ``computeMemberships`` objective.

    The selected Java loop visits only interior voxels in x/y/z order. For a
    voxel whose denominator is strictly greater than zero, it adds each
    ``class_energy[k] * membership[k] * membership[k]`` in ascending class
    order. After that conditional block, it always adds the lesion and
    black-hole energy terms, even when the denominator is zero, negative, or
    NaN. All arithmetic and every accumulator update use Java-compatible
    float32 rounding.

    Inputs are already-computed energy terms and already-normalized maps. This
    function does not calculate the energy terms, memberships, or denominator.
    It accepts IEEE NaN and infinity values and preserves their arithmetic
    behavior. The returned ``running_total`` records the cumulative Java
    accumulator after each interior voxel; zero outside that traversal region.
    """
    memberships = np.asarray(normalized_class_memberships, dtype=np.float32)
    energies = np.asarray(class_energies, dtype=np.float32)
    den = np.asarray(denominator, dtype=np.float32)
    lesion = np.asarray(lesion_membership, dtype=np.float32)
    black_hole = np.asarray(black_hole_membership, dtype=np.float32)
    lesion_e = np.asarray(lesion_energy, dtype=np.float32)
    black_hole_e = np.asarray(black_hole_energy, dtype=np.float32)

    if memberships.ndim != 4 or memberships.shape[-1] == 0:
        raise ValueError("normalized_class_memberships must have shape (x, y, z, classes)")
    spatial_shape = memberships.shape[:3]
    if min(spatial_shape) < 3:
        raise ValueError("each spatial membership dimension must contain an interior and boundary")
    if energies.shape != memberships.shape:
        raise ValueError("class_energies must match normalized_class_memberships")
    if any(value.shape != spatial_shape for value in (den, lesion, black_hole, lesion_e, black_hole_e)):
        raise ValueError("denominator and auxiliary maps must match the membership spatial shape")

    running = np.zeros(spatial_shape, dtype=np.float32)
    total = np.float32(0.0)
    with np.errstate(invalid="ignore", over="ignore", under="ignore"):
        for x in range(1, spatial_shape[0] - 1):
            for y in range(1, spatial_shape[1] - 1):
                for z in range(1, spatial_shape[2] - 1):
                    if den[x, y, z] > np.float32(0.0):
                        for k in range(memberships.shape[-1]):
                            value = np.float32(energies[x, y, z, k] * memberships[x, y, z, k])
                            value = np.float32(value * memberships[x, y, z, k])
                            total = np.float32(total + value)
                    lesion_term = np.float32(lesion_e[x, y, z] * lesion[x, y, z])
                    lesion_term = np.float32(lesion_term * lesion[x, y, z])
                    black_hole_term = np.float32(black_hole_e[x, y, z] * black_hole[x, y, z])
                    black_hole_term = np.float32(black_hole_term * black_hole[x, y, z])
                    auxiliary = np.float32(lesion_term + black_hole_term)
                    total = np.float32(total + auxiliary)
                    running[x, y, z] = total
    return LesionObjectiveReduction(total=total, running_total=running)
