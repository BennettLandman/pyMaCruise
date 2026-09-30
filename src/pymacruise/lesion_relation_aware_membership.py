"""Use the selected relation map when updating LesionTOADS memberships.

The source selects direct versus prior membership per voxel and class through
``relationList[relationMap[x][y][z]][k]``. P39zzzzs currently accepts a single
relation type per class; this adapter reuses its validated one-voxel update
kernel while selecting that type from the supplied relation state.
Paper context: Shiee et al. (2010), ``docs/papers/lesion-toads.md`` describes
modified lesion-aware relationship terms. This adapter reproduces the selected
source relation lookup for a bounded update, not the full topology evolution.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_objective import reduce_lesion_toads_objective
from .lesion_objective_membership_pass import (
    LesionObjectiveMembershipPass,
    run_lesion_toads_objective_membership_pass,
)
from .lesion_relations import LesionRelationList


@dataclass(frozen=True)
class LesionRelationAwareMembership:
    """Membership output and source-derived per-voxel/class dispatch."""

    membership: LesionObjectiveMembershipPass
    direct_dispatch: NDArray[np.bool_]


def run_lesion_toads_relation_aware_membership_pass(
    initial_class_memberships: ArrayLike,
    initial_lesion_membership: ArrayLike,
    initial_black_hole_membership: ArrayLike,
    initial_outlier_membership: ArrayLike,
    class_data_terms: ArrayLike,
    lesion_data_terms: ArrayLike,
    black_hole_data_terms: ArrayLike,
    outmask_data_terms: ArrayLike,
    shape_values: ArrayLike,
    smoothing_priors: ArrayLike,
    class_types: Sequence[str],
    relations: LesionRelationList,
    *,
    smoothing: float,
    shape_priors: ArrayLike,
    class_similarities: ArrayLike,
    lesion_similarities: ArrayLike,
    black_hole_similarities: ArrayLike,
    prior_coefficient: float,
    zero: float = 1.0e-10,
    infinity: float = 1.0e10,
) -> LesionRelationAwareMembership:
    """Run one interior update with the selected per-voxel relation lookup.

    Each source voxel is evaluated in x/y/z order. A 3-cube patch is passed to
    the already validated fixed-input membership kernel, which evaluates only
    its center; subsequent voxels therefore see the same in-place updates as
    the Java sweep. ``relations`` must be the relation-list state left by the
    continuing wrapper branch. The caller controls whether that state is
    intentionally stale relative to segmentation after subsequent thinning.
    """
    current = np.asarray(initial_class_memberships, dtype=np.float32).copy()
    lesion = np.asarray(initial_lesion_membership, dtype=np.float32).copy()
    black_hole = np.asarray(initial_black_hole_membership, dtype=np.float32).copy()
    outlier = np.asarray(initial_outlier_membership, dtype=np.float32).copy()
    class_data = np.asarray(class_data_terms, dtype=np.float32)
    lesion_data = np.asarray(lesion_data_terms, dtype=np.float32)
    black_data = np.asarray(black_hole_data_terms, dtype=np.float32)
    outmask_data = np.asarray(outmask_data_terms, dtype=np.float32)
    shapes = np.asarray(shape_values, dtype=np.float32)
    smooth_priors = np.asarray(smoothing_priors, dtype=np.float32)
    atlas_priors = np.asarray(shape_priors, dtype=np.float32)
    similarities = np.asarray(class_similarities, dtype=np.float32)
    lesion_similarity = np.asarray(lesion_similarities, dtype=np.float32)
    black_similarity = np.asarray(black_hole_similarities, dtype=np.float32)
    if not isinstance(relations, LesionRelationList):
        raise TypeError("relations must be a LesionRelationList")
    if current.ndim != 4 or min(current.shape[:3]) < 3 or current.shape[-1] == 0:
        raise ValueError("memberships must have shape (x>=3,y>=3,z>=3,classes)")
    spatial, classes = current.shape[:3], current.shape[-1]
    types = tuple(str(value).upper() for value in class_types)
    if len(types) != classes:
        raise ValueError("one class type is required per membership channel")
    if relations.relation_map.shape != spatial:
        raise ValueError("relation map must match the membership spatial grid")
    if relations.relation_membership.shape[1] != classes:
        raise ValueError("relation membership table must match the membership class axis")
    if class_data.shape != current.shape or shapes.shape != current.shape:
        raise ValueError("class data and shape maps must match memberships")
    if any(value.shape != spatial for value in (lesion, black_hole, outlier, lesion_data, black_data, outmask_data)):
        raise ValueError("all auxiliary maps must match the membership spatial grid")
    if smooth_priors.shape != (classes,) or atlas_priors.shape != (classes,):
        raise ValueError("smoothing and atlas priors must have one value per class")
    if similarities.shape != (classes, classes) or lesion_similarity.shape != (classes,) or black_similarity.shape != (classes,):
        raise ValueError("similarity inputs must match the class dimensions")
    if not np.isfinite(current).all() or not np.isfinite(lesion).all() or not np.isfinite(black_hole).all() or not np.isfinite(outlier).all():
        raise ValueError("initial membership maps must be finite")

    class_energy = np.zeros_like(current, dtype=np.float32)
    lesion_energy = np.zeros(spatial, dtype=np.float32)
    black_energy = np.zeros(spatial, dtype=np.float32)
    denominator = np.zeros(spatial, dtype=np.float32)
    direct_dispatch = np.zeros(current.shape, dtype=np.bool_)
    nx, ny, nz = spatial
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                relation_id = int(relations.relation_map[x, y, z])
                if not 0 <= relation_id < relations.relation_membership.shape[0]:
                    raise ValueError(f"relation map ID {relation_id} has no relation row")
                direct = relations.relation_membership[relation_id]
                direct_dispatch[x, y, z] = direct
                relation_types = tuple("direct" if bool(direct[k]) else "prior" for k in range(classes))
                slices = (slice(x - 1, x + 2), slice(y - 1, y + 2), slice(z - 1, z + 2))
                voxel = run_lesion_toads_objective_membership_pass(
                    current[slices],
                    lesion[slices],
                    black_hole[slices],
                    outlier[slices],
                    class_data[slices],
                    lesion_data[slices],
                    black_data[slices],
                    outmask_data[slices],
                    shapes[slices],
                    smooth_priors,
                    types,
                    relation_types,
                    smoothing=smoothing,
                    shape_priors=atlas_priors,
                    class_similarities=similarities,
                    lesion_similarities=lesion_similarity,
                    black_hole_similarities=black_similarity,
                    prior_coefficient=prior_coefficient,
                    zero=zero,
                    infinity=infinity,
                )
                center = (1, 1, 1)
                current[x, y, z] = voxel.class_memberships[center]
                lesion[x, y, z] = voxel.lesion_membership[center]
                black_hole[x, y, z] = voxel.black_hole_membership[center]
                outlier[x, y, z] = voxel.outlier_membership[center]
                denominator[x, y, z] = voxel.denominator[center]
                class_energy[x, y, z] = voxel.class_energy[center]
                lesion_energy[x, y, z] = voxel.lesion_energy[center]
                black_energy[x, y, z] = voxel.black_hole_energy[center]

    reduced = reduce_lesion_toads_objective(
        current, class_energy, denominator, lesion, black_hole, lesion_energy, black_energy
    )
    output = LesionObjectiveMembershipPass(
        class_memberships=current,
        lesion_membership=lesion,
        black_hole_membership=black_hole,
        outlier_membership=outlier,
        denominator=denominator,
        class_energy=class_energy,
        lesion_energy=lesion_energy,
        black_hole_energy=black_energy,
        objective=reduced.total,
        running_objective=reduced.running_total,
    )
    return LesionRelationAwareMembership(output, direct_dispatch)
