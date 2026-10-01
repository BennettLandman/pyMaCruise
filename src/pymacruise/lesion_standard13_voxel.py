"""One-voxel composition for the selected 13-class LesionTOADS profile.

This fixed-input slice composes selected VandyCRUISE kernels in the atlas
order Background, Sulcal-CSF, Ventricles, Cerebellum-GM, Cerebrum-GM,
Caudate, Thalamus, Putamen, Brainstem, Cerebellum-WM,
OuterVentricular-WM, InterVentricular-WM, Frontal-WM.  OPTIMIZED_CSF follows
the source's active default ``outlier`` data-energy fallback because that
dedicated branch is commented out.  The caller supplies the six face
memberships and transformed atlas state; this function does not generate
registration, topology, or distance maps, nor iterate a full volume.

The MaCRUISEplus invocation values are pinned: first atlas coefficient 2.0,
prior scale 0.1, smoothing 0.2, and distance limits 3/25/2. Constructor
scaling by 13 classes and 14 lesion-aware smoothing slots is performed here
with Java-compatible float32 steps.  The first initialization call has
``gm_flag=False`` in the selected wrapper sequence.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_data_terms import compute_lesion_toads_unfielded_data_energies
from .lesion_membership import (
    convert_lesion_toads_class_membership,
    update_lesion_toads_prior_energies,
    update_lesion_toads_smoothing_energies,
    modulate_lesion_toads_membership_by_distances,
)
from .lesion_similarity import compute_lesion_toads_similarity
from .lesion_objective import reduce_lesion_toads_objective


STANDARD13_NAMES = (
    "Background", "Sulcal-CSF", "Ventricles", "Cerebellum-GM", "Cerebrum-GM",
    "Caudate", "Thalamus", "Putamen", "Brainstem", "Cerebellum-WM",
    "OuterVentricular-WM", "InterVentricular-WM", "Frontal-WM",
)
STANDARD13_TYPES = (
    "MASK", "OPTIMIZED_CSF", "OPTIMIZED_CSF", "OBJECT", "OBJECT", "OBJECT",
    "OBJECT", "OBJECT", "OBJECT", "OBJECT", "LESION", "OBJECT", "OBJECT",
)
STANDARD13_CLASS_SMOOTHING = np.float32(
    np.float32(np.float32(np.float32(0.2) / np.float32(6.0)) / np.float32(14.0))
)
STANDARD13_PRIOR_COEFFICIENT = np.float32(np.float32(2.0) / np.float32(13.0))
STANDARD13_PRIOR_SCALE = np.float32(0.1)


@dataclass(frozen=True)
class Standard13VoxelResult:
    class_memberships: NDArray[np.float32]
    lesion_membership: np.float32
    black_hole_membership: np.float32
    outlier_membership: np.float32
    lesion_relation: np.float32
    class_energies: NDArray[np.float32]
    lesion_energy: np.float32
    black_hole_energy: np.float32
    denominator: np.float32
    objective: np.float32
    class_similarity: NDArray[np.float32]
    lesion_similarity: NDArray[np.float32]


def compute_lesion_toads_standard13_voxel(
    image_values: ArrayLike,
    class_centroids: ArrayLike,
    class_stddev: ArrayLike,
    lesion_centroid: ArrayLike,
    lesion_stddev: ArrayLike,
    channel_weights: ArrayLike,
    shape_values: ArrayLike,
    shape_priors: ArrayLike,
    smoothing_priors: ArrayLike,
    direct_relations: ArrayLike,
    class_memberships: ArrayLike,
    lesion_membership: float,
    black_hole_membership: float,
    outlier_membership: float,
    lesion_relation: float,
    six_face_class_memberships: ArrayLike,
    six_face_lesion_memberships: ArrayLike,
    six_face_black_hole_memberships: ArrayLike,
    *,
    outlier_energy: float,
    lesion_power_exponent: float,
    ventricle_distance: int,
    gray_matter_distance: int,
    brainstem_distance: int,
    gm_flag: bool = False,
    spread: float = 0.5,
) -> Standard13VoxelResult:
    """Update one interior voxel from a complete fixed pre-call state.

    ``image_values`` is exactly ``(T1_MPRAGE, FLAIR)``.  Face order is
    X+, X-, Y+, Y-, Z+, Z-.  ``direct_relations[k]`` corresponds to the
    selected ``relationList[relationMap[x][y][z]][k]`` lookup.  Input shape
    values are the transformed atlas probabilities at this voxel. Distance
    maps are supplied as outputs of the source's optional distance builders;
    their generation is outside this composition. Only an active positive
    shape/prior profile is accepted, matching the source guard.
    """
    image = np.asarray(image_values, dtype=np.float32)
    centers = np.asarray(class_centroids, dtype=np.float32)
    deviations = np.asarray(class_stddev, dtype=np.float32)
    lesion_center = np.asarray(lesion_centroid, dtype=np.float32)
    lesion_deviation = np.asarray(lesion_stddev, dtype=np.float32)
    weights = np.asarray(channel_weights, dtype=np.float32)
    shapes = np.asarray(shape_values, dtype=np.float32)
    priors = np.asarray(shape_priors, dtype=np.float32)
    smoothing_coefficients = np.asarray(smoothing_priors, dtype=np.float32)
    direct = np.asarray(direct_relations, dtype=np.bool_)
    current = np.asarray(class_memberships, dtype=np.float32).copy()
    neighbors = np.asarray(six_face_class_memberships, dtype=np.float32)
    lesion_faces = np.asarray(six_face_lesion_memberships, dtype=np.float32)
    black_hole_faces = np.asarray(six_face_black_hole_memberships, dtype=np.float32)
    lesion = np.float32(lesion_membership)
    black_hole = np.float32(black_hole_membership)
    outlier = np.float32(outlier_membership)
    lesion_rel = np.float32(lesion_relation)
    if not isinstance(gm_flag, (bool, np.bool_)):
        raise ValueError("gm_flag must be boolean")

    count = 13
    if image.shape != (2,) or centers.shape != (2, count) or deviations.shape != centers.shape:
        raise ValueError("standard profile requires two channels and 13 class centroids/deviations")
    if lesion_center.shape != (2,) or lesion_deviation.shape != (2,) or weights.shape != (2,):
        raise ValueError("lesion statistics and channel weights must have two T1+FLAIR values")
    for name, array, shape in (
        ("shape_values", shapes, (count,)), ("shape_priors", priors, (count,)),
        ("smoothing_priors", smoothing_coefficients, (count,)),
        ("direct_relations", direct, (count,)), ("class_memberships", current, (count,)),
        ("six_face_class_memberships", neighbors, (6, count)),
        ("six_face_lesion_memberships", lesion_faces, (6,)),
        ("six_face_black_hole_memberships", black_hole_faces, (6,)),
    ):
        if array.shape != shape:
            raise ValueError(f"{name} must have shape {shape}")
    numeric = (image, centers, deviations, lesion_center, lesion_deviation, weights, shapes,
               priors, smoothing_coefficients, current, neighbors, lesion_faces, black_hole_faces)
    if any(not np.isfinite(array).all() for array in numeric):
        raise ValueError("numeric state must contain only finite values")
    if np.any(shapes <= np.float32(0.0)):
        raise ValueError("all 13 classes must have positive transformed shape at this voxel")

    similarities = compute_lesion_toads_similarity(
        centers, deviations, lesion_center, lesion_deviation, weights,
        ("mask", "optimized_csf", "optimized_csf", "object", "object", "object",
         "object", "object", "object", "object", "lesion", "object", "object"),
        STANDARD13_NAMES, prior_scale=float(STANDARD13_PRIOR_SCALE), t1_channel=0,
        black_hole_similarity=np.zeros(count, dtype=np.float32),
    )
    class_energy = np.zeros(count, dtype=np.float32)
    lesion_energy = np.float32(0.0)
    black_hole_energy = np.float32(0.0)
    denominator = np.float32(0.0)

    for k, kind in enumerate(STANDARD13_TYPES):
        if shapes[k] <= np.float32(0.0):
            current[k] = np.float32(0.0)
            if kind == "LESION":
                lesion = black_hole = lesion_rel = np.float32(0.0)
            continue

        if kind == "MASK":
            raw_class = np.float32(0.0)
            for channel in range(2):
                mismatch = np.float32(0.0 if image[channel] == np.float32(0.0) else 1.0)
                raw_class = np.float32(raw_class + np.float32(weights[channel] * mismatch))
            raw_lesion = np.float32(0.0)
        elif kind == "OPTIMIZED_CSF":
            # The source's OPTIMIZED_CSF energy block is commented out; its
            # active final else assigns num=outlier.
            raw_class = np.float32(outlier_energy)
            raw_lesion = np.float32(0.0)
        else:
            energies = compute_lesion_toads_unfielded_data_energies(
                image, centers[:, k], lesion_center, np.zeros(2, dtype=np.float32), weights,
                object_type="lesion" if kind == "LESION" else "object",
                flair_channel=1, use_lesion_weight=(kind == "LESION"),
                lesion_power_exponent=(float(np.float32(lesion_power_exponent)) if kind == "LESION" else None),
                cached_lesion_power_exponent=(float(np.float32(lesion_power_exponent)) if kind == "LESION" else None),
                structure_name=STANDARD13_NAMES[k], t1_channel=0,
            )
            raw_class = energies.class_energy
            raw_lesion = energies.lesion_energy if kind == "LESION" else np.float32(0.0)

        smoothed = update_lesion_toads_smoothing_energies(
            neighbors, lesion_faces, black_hole_faces, class_index=k,
            class_type="MASK" if kind == "MASK" else ("LESION" if kind == "LESION" else "OBJECT"),
            smoothing_prior=float(smoothing_coefficients[k]),
            smoothing=float(STANDARD13_CLASS_SMOOTHING),
            class_membership_data_term=float(raw_class),
            lesion_membership_data_term=float(raw_lesion),
            black_hole_membership_data_term=0.0, outmask_membership_data_term=0.0,
            class_energy=float(raw_class), lesion_energy=float(raw_lesion), black_hole_energy=0.0,
        )
        cls_term = smoothed.class_membership_data_term
        lesion_term = smoothed.lesion_membership_data_term
        cls_energy = smoothed.class_energy
        lesion_energy_value = smoothed.lesion_energy
        prior_update = update_lesion_toads_prior_energies(
            priors, similarities.class_similarity[k], similarities.lesion_similarity,
            similarities.black_hole_similarity, shapes, class_index=k,
            class_type="MASK" if kind == "MASK" else ("LESION" if kind == "LESION" else "OBJECT"),
            prior_coefficient=float(STANDARD13_PRIOR_COEFFICIENT),
            class_membership_data_term=float(cls_term),
            lesion_membership_data_term=float(lesion_term), black_hole_membership_data_term=0.0,
            outmask_membership_data_term=0.0, class_energy=float(cls_energy),
            lesion_energy=float(lesion_energy_value), black_hole_energy=0.0,
            has_lesion=True, lesion_class_index=10, lesion_class_type="LESION",
        )
        relation = "direct" if direct[k] else "prior"
        current[k] = convert_lesion_toads_class_membership(
            float(prior_update.class_membership_data_term), float(shapes[k]), relation_type=relation,
        )
        class_energy[k] = prior_update.class_energy
        if kind == "LESION":
            lesion = convert_lesion_toads_class_membership(
                float(prior_update.lesion_membership_data_term), float(shapes[k]), relation_type=relation,
            )
            lesion_rel = np.float32(1.0) if direct[k] else np.float32(
                np.float32(2.0) * shapes[k] if np.float32(2.0) * shapes[k] < np.float32(1.0) else 1.0
            )
            lesion_energy = prior_update.lesion_energy
            modulated = modulate_lesion_toads_membership_by_distances(
                float(lesion), float(lesion_rel), ventricle_distance=ventricle_distance,
                gray_matter_distance=gray_matter_distance, brainstem_distance=brainstem_distance,
    max_ventricle_distance=2, max_gray_matter_distance=3, max_brainstem_distance=25,
                gm_flag=gm_flag, spread=float(spread),
            )
            lesion, lesion_rel = modulated.membership, modulated.relation
        elif STANDARD13_NAMES[k] == "Cerebrum-GM" and ventricle_distance != -1:
            current[k] = np.float32(current[k] * np.float32(ventricle_distance) / np.float32(2.0))

        if kind == "LESION":
            denominator = np.float32(denominator + lesion)
        denominator = np.float32(denominator + current[k])

    if denominator > np.float32(0.0):
        for k in range(count):
            current[k] = np.float32(current[k] / denominator)
        lesion = np.float32(lesion / denominator)
        # Source's class-normalization loop leaves BH/outlier unchanged because
        # the selected 13-class atlas has neither BLACKHOLE nor OUTLIER classes.

    # Reuse the source-ordered objective reducer over a one-voxel interior.
    state = np.zeros((3, 3, 3, count), dtype=np.float32)
    energies = np.zeros_like(state)
    lesion_map = np.zeros((3, 3, 3), dtype=np.float32)
    bh_map = np.zeros_like(lesion_map)
    denominator_map = np.zeros_like(lesion_map)
    lesion_energy_map = np.zeros_like(lesion_map)
    bh_energy_map = np.zeros_like(lesion_map)
    state[1, 1, 1] = current
    energies[1, 1, 1] = class_energy
    lesion_map[1, 1, 1] = lesion
    bh_map[1, 1, 1] = black_hole
    denominator_map[1, 1, 1] = denominator
    lesion_energy_map[1, 1, 1] = lesion_energy
    bh_energy_map[1, 1, 1] = black_hole_energy
    objective = reduce_lesion_toads_objective(
        state, energies, denominator_map, lesion_map, bh_map, lesion_energy_map, bh_energy_map,
    ).total
    return Standard13VoxelResult(
        class_memberships=current, lesion_membership=lesion,
        black_hole_membership=black_hole, outlier_membership=outlier,
        lesion_relation=lesion_rel, class_energies=class_energy,
        lesion_energy=lesion_energy, black_hole_energy=black_hole_energy,
        denominator=denominator, objective=objective,
        class_similarity=similarities.class_similarity,
        lesion_similarity=similarities.lesion_similarity,
    )
