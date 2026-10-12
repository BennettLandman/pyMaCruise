"""Numba implementation of the selected source-ordered standard-13 sweep.

The kernel keeps the Java x/y/z/class traversal and float32 update order while
moving the hot scalar loops out of Python. Input validation and centroid
similarity construction remain in :mod:`lesion_standard13_volume`.
"""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True, fastmath=False)
def standard13_membership_sweep_numba(
    images: np.ndarray,
    class_centroids: np.ndarray,
    channel_weights: np.ndarray,
    lesion_centroid: np.ndarray,
    lesion_energy_map: np.ndarray,
    shape_values: np.ndarray,
    shape_priors: np.ndarray,
    smoothing_priors: np.ndarray,
    direct_relations: np.ndarray,
    class_similarity: np.ndarray,
    lesion_similarity: np.ndarray,
    initial_class_memberships: np.ndarray,
    initial_lesion_membership: np.ndarray,
    initial_black_hole_membership: np.ndarray,
    initial_outlier_membership: np.ndarray,
    initial_lesion_relation: np.ndarray,
    ventricle_distances: np.ndarray,
    gray_matter_distances: np.ndarray,
    brainstem_distances: np.ndarray,
    outlier_energy: np.float32,
    gm_flag: bool,
    spread: np.float32,
    source_zero_variance_nan: bool,
    class_smoothing: np.float32,
    prior_coefficient: np.float32,
) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.float32, np.ndarray,
]:
    """Apply one exact ordered volume update for the fixed 13-class atlas."""
    current = initial_class_memberships.copy()
    lesion = initial_lesion_membership.copy()
    black_hole = initial_black_hole_membership.copy()
    outlier = initial_outlier_membership.copy()
    lesion_relation = initial_lesion_relation.copy()
    spatial = current.shape[:3]
    class_energy = np.zeros(current.shape, dtype=np.float32)
    lesion_energy = np.zeros(spatial, dtype=np.float32)
    black_hole_energy = np.zeros(spatial, dtype=np.float32)
    denominator = np.zeros(spatial, dtype=np.float32)
    class_energy_state = np.zeros(13, dtype=np.float32)
    lesion_energy_state = np.float32(0.0)

    for x in range(1, spatial[0] - 1):
        for y in range(1, spatial[1] - 1):
            for z in range(1, spatial[2] - 1):
                # Snapshot the six faces once for this voxel. Each following
                # class update reads the same pre-class-loop neighborhood.
                faces = np.empty((6, 13), dtype=np.float32)
                lesion_faces = np.empty(6, dtype=np.float32)
                black_hole_faces = np.empty(6, dtype=np.float32)
                for k in range(13):
                    faces[0, k] = current[x + 1, y, z, k]
                    faces[1, k] = current[x - 1, y, z, k]
                    faces[2, k] = current[x, y + 1, z, k]
                    faces[3, k] = current[x, y - 1, z, k]
                    faces[4, k] = current[x, y, z + 1, k]
                    faces[5, k] = current[x, y, z - 1, k]
                lesion_faces[0] = lesion[x + 1, y, z]
                lesion_faces[1] = lesion[x - 1, y, z]
                lesion_faces[2] = lesion[x, y + 1, z]
                lesion_faces[3] = lesion[x, y - 1, z]
                lesion_faces[4] = lesion[x, y, z + 1]
                lesion_faces[5] = lesion[x, y, z - 1]
                black_hole_faces[0] = black_hole[x + 1, y, z]
                black_hole_faces[1] = black_hole[x - 1, y, z]
                black_hole_faces[2] = black_hole[x, y + 1, z]
                black_hole_faces[3] = black_hole[x, y - 1, z]
                black_hole_faces[4] = black_hole[x, y, z + 1]
                black_hole_faces[5] = black_hole[x, y, z - 1]

                voxel_shapes = shape_values[x, y, z]
                denominator_value = np.float32(0.0)
                for k in range(13):
                    kind = 0 if k == 0 else (1 if k == 1 or k == 2 else (2 if k == 10 else 3))
                    shape = voxel_shapes[k]
                    if not shape > np.float32(0.0):
                        current[x, y, z, k] = np.float32(0.0)
                        continue

                    if kind == 0:
                        raw_class = np.float32(0.0)
                        for channel in range(2):
                            if images[x, y, z, channel] != np.float32(0.0):
                                raw_class = np.float32(raw_class + channel_weights[channel])
                        raw_lesion = np.float32(0.0)
                    elif kind == 1:
                        raw_class = outlier_energy
                        raw_lesion = np.float32(0.0)
                    else:
                        raw_class = np.float32(0.0)
                        for channel in range(2):
                            difference = np.float32(
                                images[x, y, z, channel] - class_centroids[channel, k]
                            )
                            weighted_difference = np.float32(channel_weights[channel] * difference)
                            contribution = np.float32(weighted_difference * difference)
                            raw_class = np.float32(raw_class + contribution)
                        raw_lesion = lesion_energy_map[x, y, z] if kind == 2 else np.float32(0.0)

                    other_neighbor_sum = np.float32(0.0)
                    lesion_neighbor_sum = np.float32(0.0)
                    black_hole_neighbor_sum = np.float32(0.0)
                    for face in range(6):
                        for class_no in range(13):
                            if class_no != k:
                                value = faces[face, class_no]
                                square = np.float32(value * value)
                                other_neighbor_sum = np.float32(other_neighbor_sum + square)
                        value = lesion_faces[face]
                        lesion_neighbor_sum = np.float32(
                            lesion_neighbor_sum + np.float32(value * value)
                        )
                        value = black_hole_faces[face]
                        black_hole_neighbor_sum = np.float32(
                            black_hole_neighbor_sum + np.float32(value * value)
                        )

                    current_neighbor_sum = np.float32(0.0)
                    if kind == 2:
                        for face in (1, 0, 3, 2, 5, 4):
                            value = faces[face, k]
                            current_neighbor_sum = np.float32(
                                current_neighbor_sum + np.float32(value * value)
                            )

                    smoothing_prior = smoothing_priors[k]
                    all_neighbor_sum = np.float32(other_neighbor_sum + lesion_neighbor_sum)
                    all_neighbor_sum = np.float32(all_neighbor_sum + black_hole_neighbor_sum)
                    class_full = np.float32(
                        np.float32(smoothing_prior * class_smoothing) * all_neighbor_sum
                    )
                    half_prior = np.float32(np.float32(0.5) * smoothing_prior)
                    class_half = np.float32(
                        np.float32(half_prior * class_smoothing) * all_neighbor_sum
                    )
                    class_term = np.float32(raw_class + class_full)
                    class_energy_value = np.float32(raw_class + class_half)
                    lesion_term = raw_lesion
                    lesion_energy_value = lesion_energy_state
                    if kind == 2:
                        lesion_sum_for_update = np.float32(other_neighbor_sum + current_neighbor_sum)
                        lesion_sum_for_update = np.float32(
                            lesion_sum_for_update + black_hole_neighbor_sum
                        )
                        lesion_full = np.float32(
                            np.float32(smoothing_prior * class_smoothing)
                            * lesion_sum_for_update
                        )
                        lesion_half = np.float32(
                            np.float32(half_prior * class_smoothing) * lesion_sum_for_update
                        )
                        lesion_term = np.float32(lesion_term + lesion_full)
                        lesion_energy_value = np.float32(raw_lesion + lesion_half)

                    class_prior = np.float32(0.0)
                    lesion_prior = np.float32(0.0)
                    for m in range(13):
                        other_shape = voxel_shapes[m]
                        if m != k:
                            product = np.float32(shape_priors[m] * class_similarity[k, m])
                            product = np.float32(product * other_shape)
                            product = np.float32(product * other_shape)
                            class_prior = np.float32(class_prior + product)
                        if kind == 2:
                            product = np.float32(shape_priors[m] * lesion_similarity[m])
                            product = np.float32(product * other_shape)
                            product = np.float32(product * other_shape)
                            lesion_prior = np.float32(lesion_prior + product)
                    lesion_shape = voxel_shapes[10]
                    product = np.float32(shape_priors[10] * lesion_similarity[k])
                    product = np.float32(product * lesion_shape)
                    product = np.float32(product * lesion_shape)
                    class_prior = np.float32(class_prior + product)

                    class_term = np.float32(
                        class_term + np.float32(prior_coefficient * class_prior)
                    )
                    prior_half = np.float32(
                        np.float32(np.float32(0.5) * prior_coefficient) * class_prior
                    )
                    class_energy_value = np.float32(class_energy_value + prior_half)
                    if kind == 2:
                        lesion_term = np.float32(
                            lesion_term + np.float32(prior_coefficient * lesion_prior)
                        )
                        lesion_half = np.float32(
                            np.float32(np.float32(0.5) * prior_coefficient) * lesion_prior
                        )
                        lesion_energy_value = np.float32(lesion_energy_value + lesion_half)

                    relation_is_direct = direct_relations[x, y, z, k]
                    energy_value = class_term
                    if source_zero_variance_nan:
                        reciprocal = (
                            np.float32(1.0) / energy_value
                            if energy_value > np.float32(1.0e-10)
                            else np.float32(1.0e10)
                        )
                    else:
                        reciprocal = (
                            np.float32(1.0) / energy_value
                            if energy_value > np.float32(1.0e-10)
                            else np.float32(1.0e10)
                        )
                    if not relation_is_direct:
                        doubled_shape = np.float32(np.float32(2.0) * shape)
                        factor = np.float32(1.0) if np.float32(1.0) < doubled_shape else doubled_shape
                        reciprocal = np.float32(factor * reciprocal)
                    current[x, y, z, k] = reciprocal
                    class_energy_state[k] = class_energy_value
                    class_energy[x, y, z, k] = class_energy_value

                    if kind == 2:
                        relation = np.float32(1.0) if relation_is_direct else np.float32(
                            np.float32(2.0) * shape
                            if np.float32(2.0) * shape < np.float32(1.0)
                            else 1.0
                        )
                        lesion_value = (
                            np.float32(1.0) / lesion_term
                            if lesion_term > np.float32(1.0e-10)
                            else np.float32(1.0e10)
                        )
                        if not relation_is_direct:
                            doubled_shape = np.float32(np.float32(2.0) * shape)
                            factor = np.float32(1.0) if np.float32(1.0) < doubled_shape else doubled_shape
                            lesion_value = np.float32(factor * lesion_value)

                        vent = ventricle_distances[x, y, z]
                        gray = gray_matter_distances[x, y, z]
                        brain = brainstem_distances[x, y, z]
                        if vent != -1:
                            numerator = spread if vent == 0 else np.float32(vent)
                            lesion_value = np.float32(
                                np.float32(numerator * lesion_value) / np.float32(2.0)
                            )
                            factor = (
                                np.float32(numerator / np.float32(2.0))
                                if vent == 0
                                else np.float32(vent // 2)
                            )
                            relation = np.float32(relation * factor)
                        elif gm_flag and gray != -1:
                            numerator = spread if gray == 0 else np.float32(gray)
                            lesion_value = np.float32(
                                np.float32(numerator * lesion_value) / np.float32(3.0)
                            )
                            factor = (
                                np.float32(numerator / np.float32(3.0))
                                if gray == 0
                                else np.float32(gray // 3)
                            )
                            relation = np.float32(relation * factor)
                        if brain != -1:
                            lesion_value = np.float32(
                                np.float32(np.float32(brain) * lesion_value) / np.float32(25.0)
                            )
                            relation = np.float32(relation * np.float32(brain // 25))
                        lesion[x, y, z] = lesion_value
                        lesion_relation[x, y, z] = relation
                        lesion_energy_state = lesion_energy_value
                        lesion_energy[x, y, z] = lesion_energy_value
                    elif k == 4:
                        vent = ventricle_distances[x, y, z]
                        if vent != -1:
                            current[x, y, z, k] = np.float32(
                                np.float32(current[x, y, z, k] * np.float32(vent))
                                / np.float32(2.0)
                            )

                    if kind == 2:
                        denominator_value = np.float32(denominator_value + lesion[x, y, z])
                    denominator_value = np.float32(
                        denominator_value + current[x, y, z, k]
                    )

                for k in range(13):
                    class_energy[x, y, z, k] = class_energy_state[k]
                lesion_energy[x, y, z] = lesion_energy_state
                denominator[x, y, z] = denominator_value
                if denominator_value > np.float32(0.0):
                    for k in range(13):
                        current[x, y, z, k] = np.float32(
                            current[x, y, z, k] / denominator_value
                        )
                    lesion[x, y, z] = np.float32(lesion[x, y, z] / denominator_value)
                black_hole_energy[x, y, z] = np.float32(0.0)

    running = np.zeros(spatial, dtype=np.float32)
    total = np.float32(0.0)
    for x in range(1, spatial[0] - 1):
        for y in range(1, spatial[1] - 1):
            for z in range(1, spatial[2] - 1):
                if denominator[x, y, z] > np.float32(0.0):
                    for k in range(13):
                        value = np.float32(class_energy[x, y, z, k] * current[x, y, z, k])
                        value = np.float32(value * current[x, y, z, k])
                        total = np.float32(total + value)
                lesion_term = np.float32(lesion_energy[x, y, z] * lesion[x, y, z])
                lesion_term = np.float32(lesion_term * lesion[x, y, z])
                black_hole_term = np.float32(
                    black_hole_energy[x, y, z] * black_hole[x, y, z]
                )
                black_hole_term = np.float32(black_hole_term * black_hole[x, y, z])
                total = np.float32(total + np.float32(lesion_term + black_hole_term))
                running[x, y, z] = total

    return (
        current, lesion, black_hole, outlier, lesion_relation, denominator,
        class_energy, lesion_energy, black_hole_energy, total, running,
    )
