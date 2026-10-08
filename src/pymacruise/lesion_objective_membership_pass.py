"""Objective-reporting LesionTOADS membership pass for supported fixed inputs.

This composes the selected LesionTOADS interior update for OBJECT, LESION,
BLACKHOLE, OUTMASK, MASK, and OUTLIER class families.
Images and statistical estimation are upstream: inputs here are fixed energy,
shape, similarity, and membership maps. Every interior class must be active
(positive shape) in the safe default mode. ``source_faithful=True`` also
preserves the selected Java scan's zero-shape skip behavior: active classes
reset their energy from the current data term, while a skipped class retains
the last value in the method-local energy array for the objective.
Paper context: Shiee et al. (2010), ``docs/papers/lesion-toads.md``. This
source-selected membership/objective slice uses fixed maps and excludes
statistics estimation and unsupported skipped-shape cases.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_membership import (
    convert_lesion_toads_class_membership,
    iter_lesion_toads_six_face_neighborhoods,
    update_lesion_toads_outmask_membership,
    update_lesion_toads_prior_energies,
    update_lesion_toads_smoothing_energies,
)
from .lesion_objective import reduce_lesion_toads_objective


@dataclass(frozen=True)
class LesionObjectiveMembershipPass:
    """Updated RAS arrays and selected objective outputs."""

    class_memberships: NDArray[np.float32]
    lesion_membership: NDArray[np.float32]
    black_hole_membership: NDArray[np.float32]
    outlier_membership: NDArray[np.float32]
    denominator: NDArray[np.float32]
    class_energy: NDArray[np.float32]
    lesion_energy: NDArray[np.float32]
    black_hole_energy: NDArray[np.float32]
    objective: np.float32
    running_objective: NDArray[np.float32]


def run_lesion_toads_objective_membership_pass(
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
    relation_types: Sequence[str],
    *,
    smoothing: float,
    shape_priors: ArrayLike,
    class_similarities: ArrayLike,
    lesion_similarities: ArrayLike,
    black_hole_similarities: ArrayLike,
    prior_coefficient: float,
    zero: float = 1.0e-10,
    infinity: float = 1.0e10,
    source_faithful: bool = False,
) -> LesionObjectiveMembershipPass:
    """Run one fixed-input, source-ordered OBJECT/LESION/BH/OUTMASK pass.

    Arrays use ``(x, y, z[, classes])`` RAS index order. The function visits
    interior x/y/z coordinates and ascending classes, mutating the working
    memberships in place so later voxels observe earlier updates. Exactly one
    BLACKHOLE and OUTMASK class and at most one LESION class are supported;
    only direct and prior relations are supported. By default, all classes
    must have positive shape at every interior voxel. With
    ``source_faithful=True``, zero-shape classes are skipped as in the source;
    their membership is set to zero while their energy value retains its
    previous value. Active class and auxiliary energies start from the current
    voxel's data terms and receive same-voxel smoothing and prior additions.

    ``class_data_terms`` contains each class's raw ``num`` term,
    ``lesion_data_terms``/``black_hole_data_terms`` contain the auxiliary
    raw terms, and ``outmask_data_terms`` contains raw ``numO``. The returned
    energy maps include the selected half-weighted smoothing/prior terms.
    """
    current = np.asarray(initial_class_memberships, dtype=np.float32).copy()
    lesion = np.asarray(initial_lesion_membership, dtype=np.float32).copy()
    black_hole = np.asarray(initial_black_hole_membership, dtype=np.float32).copy()
    outlier = np.asarray(initial_outlier_membership, dtype=np.float32).copy()
    class_data = np.asarray(class_data_terms, dtype=np.float32)
    lesion_data = np.asarray(lesion_data_terms, dtype=np.float32)
    black_hole_data = np.asarray(black_hole_data_terms, dtype=np.float32)
    outmask_data = np.asarray(outmask_data_terms, dtype=np.float32)
    shapes = np.asarray(shape_values, dtype=np.float32)
    smooth_prior = np.asarray(smoothing_priors, dtype=np.float32)
    atlas_prior = np.asarray(shape_priors, dtype=np.float32)
    similarities = np.asarray(class_similarities, dtype=np.float32)
    lesion_similarity = np.asarray(lesion_similarities, dtype=np.float32)
    black_hole_similarity = np.asarray(black_hole_similarities, dtype=np.float32)

    if current.ndim != 4 or current.shape[-1] < 1 or min(current.shape[:3]) < 3:
        raise ValueError("initial_class_memberships must have shape (x>=3,y>=3,z>=3,classes>=1)")
    spatial = current.shape[:3]
    classes = current.shape[-1]
    types = tuple(class_types)
    relations = tuple(relation_types)
    if len(types) != classes or types.count("BLACKHOLE") > 1 or types.count("OUTMASK") > 1:
        raise ValueError("class_types must match classes and contain at most one BLACKHOLE and OUTMASK")
    if sum(kind == "LESION" for kind in types) > 1 or any(
        kind not in {"OBJECT", "LESION", "BLACKHOLE", "OUTMASK", "MASK", "OUTLIER"} for kind in types
    ):
        raise ValueError("unsupported class family in objective membership pass")
    if not any(kind in {"LESION", "BLACKHOLE"} for kind in types):
        raise ValueError("the objective membership pass requires a lesion-related class")
    if len(relations) != classes or any(kind not in {"direct", "prior"} for kind in relations):
        raise ValueError("relation_types must contain one direct/prior relation per class")
    if class_data.shape != current.shape or shapes.shape != current.shape:
        raise ValueError("class_data_terms and shape_values must match membership shape")
    if any(array.shape != spatial for array in (lesion, black_hole, outlier, lesion_data, black_hole_data, outmask_data)):
        raise ValueError("auxiliary input maps must match the spatial grid")
    if smooth_prior.shape != (classes,) or atlas_prior.shape != (classes,):
        raise ValueError("smoothing_priors and shape_priors must have one value per class")
    if similarities.shape != (classes, classes) or lesion_similarity.shape != (classes,) or black_hole_similarity.shape != (classes,):
        raise ValueError("similarity inputs must match the class dimensions")
    finite_arrays = (
        current, lesion, black_hole, outlier, class_data, lesion_data,
        black_hole_data, outmask_data, shapes, smooth_prior, atlas_prior,
        similarities, lesion_similarity, black_hole_similarity,
    )
    if any(not np.isfinite(array).all() for array in finite_arrays):
        raise ValueError("all map, prior, and similarity inputs must be finite")
    if not isinstance(source_faithful, (bool, np.bool_)):
        raise ValueError("source_faithful must be boolean")
    if np.any(smooth_prior < 0) or np.any(atlas_prior < 0):
        raise ValueError("smoothing_priors and shape_priors must be nonnegative")
    smooth = np.float32(smoothing)
    coefficient = np.float32(prior_coefficient)
    zero_value = np.float32(zero)
    infinity_value = np.float32(infinity)
    if not np.isfinite(smooth) or smooth < 0 or not np.isfinite(coefficient) or coefficient < 0:
        raise ValueError("smoothing and prior_coefficient must be finite and nonnegative")
    if not np.isfinite(zero_value) or zero_value < 0 or not np.isfinite(infinity_value) or infinity_value < 0:
        raise ValueError("zero and infinity must be finite and nonnegative")
    interior_shapes = shapes[1:-1, 1:-1, 1:-1]
    if np.any(interior_shapes < np.float32(0.0)) or (
        not source_faithful and np.any(interior_shapes == np.float32(0.0))
    ):
        raise ValueError("all interior classes must be active with strictly positive shape values")

    lesion_class = max((k for k, kind in enumerate(types) if kind in {"LESION", "BLACKHOLE"}), default=-1)
    if lesion_class < 0:
        raise ValueError("the supported objective pass requires a lesion-related class")
    class_energy = np.zeros(current.shape, dtype=np.float32)
    lesion_energy = np.zeros(spatial, dtype=np.float32)
    black_hole_energy = np.zeros(spatial, dtype=np.float32)
    denominator_map = np.zeros(spatial, dtype=np.float32)
    class_energy_state = np.zeros(classes, dtype=np.float32)
    lesion_energy_state = np.float32(0.0)
    black_hole_energy_state = np.float32(0.0)

    for neighborhood in iter_lesion_toads_six_face_neighborhoods(current, lesion, black_hole):
        x, y, z = neighborhood.coordinate
        denominator = np.float32(0.0)
        for k, class_type in enumerate(types):
            if source_faithful and shapes[x, y, z, k] == np.float32(0.0):
                # The source's ``if (shape[k] > 0)`` skips the update but its
                # following ``else`` zeros this class membership. energy[k]
                # is not assigned on this skipped branch.
                current[x, y, z, k] = np.float32(0.0)
                class_energy[x, y, z, k] = class_energy_state[k]
                continue
            raw_class = np.float32(class_data[x, y, z, k])
            raw_lesion = np.float32(lesion_data[x, y, z]) if class_type in {"LESION", "BLACKHOLE"} else np.float32(0.0)
            raw_black_hole = np.float32(black_hole_data[x, y, z]) if class_type == "BLACKHOLE" else np.float32(0.0)
            raw_outmask = np.float32(outmask_data[x, y, z]) if class_type == "OUTMASK" else np.float32(0.0)

            # Java stores the mixed OUTMASK data energy before spatial terms;
            # the membership update separately inverts and mixes its terms.
            energy_base = raw_class
            if class_type == "OUTMASK":
                energy_base = np.float32(
                    np.float32(np.float32(1.0) - outlier[x, y, z]) * raw_class
                    + np.float32(outlier[x, y, z] * raw_outmask)
                )

            smoothing_update = update_lesion_toads_smoothing_energies(
                neighborhood.class_memberships,
                neighborhood.lesion_memberships,
                neighborhood.black_hole_memberships,
                class_index=k,
                class_type=class_type,
                smoothing_prior=float(smooth_prior[k]),
                smoothing=float(smooth),
                class_membership_data_term=float(raw_class),
                lesion_membership_data_term=float(raw_lesion),
                black_hole_membership_data_term=float(raw_black_hole),
                outmask_membership_data_term=float(raw_outmask),
                class_energy=float(energy_base),
                lesion_energy=float(raw_lesion),
                black_hole_energy=float(raw_black_hole),
            )
            class_term = np.float32(smoothing_update.class_membership_data_term)
            lesion_term = np.float32(smoothing_update.lesion_membership_data_term)
            black_hole_term = np.float32(smoothing_update.black_hole_membership_data_term)
            outmask_term = np.float32(smoothing_update.outmask_membership_data_term)
            class_energy_value = np.float32(smoothing_update.class_energy)
            lesion_energy_value = np.float32(smoothing_update.lesion_energy)
            black_hole_energy_value = np.float32(smoothing_update.black_hole_energy)

            if coefficient > np.float32(0.0):
                prior_update = update_lesion_toads_prior_energies(
                    atlas_prior,
                    similarities[k],
                    lesion_similarity,
                    black_hole_similarity,
                    shapes[x, y, z],
                    class_index=k,
                    class_type=class_type,
                    prior_coefficient=float(coefficient),
                    class_membership_data_term=float(class_term),
                    lesion_membership_data_term=float(lesion_term),
                    black_hole_membership_data_term=float(black_hole_term),
                    outmask_membership_data_term=float(outmask_term),
                    class_energy=float(class_energy_value),
                    lesion_energy=float(lesion_energy_value),
                    black_hole_energy=float(black_hole_energy_value),
                    has_lesion=True,
                    lesion_class_index=lesion_class,
                    lesion_class_type=types[lesion_class],
                )
                class_term = np.float32(prior_update.class_membership_data_term)
                lesion_term = np.float32(prior_update.lesion_membership_data_term)
                black_hole_term = np.float32(prior_update.black_hole_membership_data_term)
                outmask_term = np.float32(prior_update.outmask_membership_data_term)
                class_energy_value = np.float32(prior_update.class_energy)
                lesion_energy_value = np.float32(prior_update.lesion_energy)
                black_hole_energy_value = np.float32(prior_update.black_hole_energy)

            membership = convert_lesion_toads_class_membership(
                float(class_term), float(shapes[x, y, z, k]),
                relation_type=relations[k], zero=float(zero_value), infinity=float(infinity_value),
            )
            if class_type in {"LESION", "BLACKHOLE"}:
                lesion_membership = convert_lesion_toads_class_membership(
                    float(lesion_term), float(shapes[x, y, z, k]),
                    relation_type=relations[k], zero=float(zero_value), infinity=float(infinity_value),
                )
                lesion[x, y, z] = lesion_membership
            if class_type == "BLACKHOLE":
                black_hole[x, y, z] = convert_lesion_toads_class_membership(
                    float(black_hole_term), float(shapes[x, y, z, k]),
                    relation_type=relations[k], zero=float(zero_value), infinity=float(infinity_value),
                )
            if class_type == "OUTMASK":
                outmask_membership = convert_lesion_toads_class_membership(
                    float(outmask_term), 1.0, relation_type="direct",
                    zero=float(zero_value), infinity=float(infinity_value),
                )
                mixed = update_lesion_toads_outmask_membership(
                    float(convert_lesion_toads_class_membership(
                        float(class_term), 1.0, relation_type="direct",
                        zero=float(zero_value), infinity=float(infinity_value),
                    )),
                    float(outmask_membership),
                    float(outlier[x, y, z]),
                    relation_type=relations[k],
                    shape_prior=float(shapes[x, y, z, k]),
                )
                membership = mixed.class_membership
                outlier[x, y, z] = mixed.outlier_membership
            current[x, y, z, k] = membership
            if class_type == "LESION":
                denominator = np.float32(denominator + lesion[x, y, z])
            elif class_type == "BLACKHOLE":
                auxiliary_sum = np.float32(lesion[x, y, z] + black_hole[x, y, z])
                denominator = np.float32(denominator + auxiliary_sum)
            denominator = np.float32(denominator + membership)

            if source_faithful:
                class_energy_state[k] = class_energy_value
                class_energy[x, y, z, k] = class_energy_state[k]
            else:
                class_energy[x, y, z, k] = class_energy_value
            if class_type in {"LESION", "BLACKHOLE"}:
                if source_faithful:
                    lesion_energy_state = lesion_energy_value
                    lesion_energy[x, y, z] = lesion_energy_state
                else:
                    lesion_energy[x, y, z] = lesion_energy_value
            if class_type == "BLACKHOLE":
                if source_faithful:
                    black_hole_energy_state = black_hole_energy_value
                    black_hole_energy[x, y, z] = black_hole_energy_state
                else:
                    black_hole_energy[x, y, z] = black_hole_energy_value

        denominator_map[x, y, z] = denominator
        if source_faithful:
            # The Java method uses its latest auxiliary energy values at this
            # voxel, including a retained value after a skipped class.
            lesion_energy[x, y, z] = lesion_energy_state
            black_hole_energy[x, y, z] = black_hole_energy_state
        if denominator > np.float32(0.0):
            for k, class_type in enumerate(types):
                current[x, y, z, k] = np.float32(current[x, y, z, k] / denominator)
                if class_type in {"LESION", "BLACKHOLE"}:
                    lesion[x, y, z] = np.float32(lesion[x, y, z] / denominator)
                if class_type == "BLACKHOLE":
                    black_hole[x, y, z] = np.float32(black_hole[x, y, z] / denominator)
                elif class_type == "OUTMASK":
                    outlier[x, y, z] = np.float32(outlier[x, y, z] / denominator)

    objective_result = reduce_lesion_toads_objective(
        current, class_energy, denominator_map, lesion, black_hole, lesion_energy, black_hole_energy
    )
    return LesionObjectiveMembershipPass(
        class_memberships=current,
        lesion_membership=lesion,
        black_hole_membership=black_hole,
        outlier_membership=outlier,
        denominator=denominator_map,
        class_energy=class_energy,
        lesion_energy=lesion_energy,
        black_hole_energy=black_hole_energy,
        objective=objective_result.total,
        running_objective=objective_result.running_total,
    )
