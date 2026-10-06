"""Volume-wide source-ordered pass for the fixed standard LesionTOADS atlas.

This composes :mod:`lesion_standard13_voxel` over interior coordinates in
Java's x/y/z order. Updated memberships are written immediately, so later
voxels read those changes through their six-face neighborhoods. Transformed
shape, relation, and optional distance state are fixed caller inputs; this
module does not construct them or claim complete MaCRUISE/LesionTOADS
equivalence.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_objective import reduce_lesion_toads_objective
from .lesion_similarity import compute_lesion_toads_similarity
from .lesion_standard13_voxel import (
    STANDARD13_NAMES,
    STANDARD13_PRIOR_SCALE,
    STANDARD13_TYPES,
    compute_lesion_toads_standard13_voxel,
)


@dataclass(frozen=True)
class Standard13VolumeResult:
    class_memberships: NDArray[np.float32]
    lesion_membership: NDArray[np.float32]
    black_hole_membership: NDArray[np.float32]
    outlier_membership: NDArray[np.float32]
    lesion_relation: NDArray[np.float32]
    denominator: NDArray[np.float32]
    class_energies: NDArray[np.float32]
    lesion_energy: NDArray[np.float32]
    black_hole_energy: NDArray[np.float32]
    objective: np.float32
    running_objective: NDArray[np.float32]


def run_lesion_toads_standard13_volume(
    images: ArrayLike,
    class_centroids: ArrayLike,
    class_stddev: ArrayLike,
    lesion_centroid: ArrayLike,
    lesion_stddev: ArrayLike,
    channel_weights: ArrayLike,
    shape_values: ArrayLike,
    shape_priors: ArrayLike,
    smoothing_priors: ArrayLike,
    direct_relations: ArrayLike,
    initial_class_memberships: ArrayLike,
    initial_lesion_membership: ArrayLike,
    initial_black_hole_membership: ArrayLike,
    initial_outlier_membership: ArrayLike,
    initial_lesion_relation: ArrayLike,
    ventricle_distances: ArrayLike,
    gray_matter_distances: ArrayLike,
    brainstem_distances: ArrayLike,
    *,
    outlier_energy: float,
    lesion_power_exponent: float,
    gm_flag: bool = False,
    spread: float = 0.5,
    source_zero_variance_nan: bool = False,
    intensity_fields: ArrayLike | None = None,
) -> Standard13VolumeResult:
    """Run the selected 13-class membership update over one volume.

    Arrays use ``(x, y, z[, classes])`` order and are assumed to share the
    supplied RAS grid. The class relation flags, transformed atlas shapes,
    and source-computed distance maps are fixed. Border voxels are untouched.
    The returned objective is reduced in the same x/y/z/class accumulation
    order as Java after all in-place voxel updates are complete.
    """
    data = np.asarray(images, dtype=np.float32)
    if data.ndim != 4 or data.shape[-1] != 2 or min(data.shape[:3], default=0) < 3:
        raise ValueError("images must have shape (x>=3,y>=3,z>=3,2) in T1_MPRAGE/FLAIR order")
    spatial = data.shape[:3]
    fields = None if intensity_fields is None else np.asarray(intensity_fields, dtype=np.float32)
    if fields is not None and (fields.shape != data.shape or not np.isfinite(fields).all()):
        raise ValueError("intensity_fields must be finite and match the image grid and channels")
    shapes = np.asarray(shape_values, dtype=np.float32)
    direct = np.asarray(direct_relations, dtype=np.bool_)
    current = np.asarray(initial_class_memberships, dtype=np.float32).copy()
    lesion = np.asarray(initial_lesion_membership, dtype=np.float32).copy()
    black_hole = np.asarray(initial_black_hole_membership, dtype=np.float32).copy()
    outlier = np.asarray(initial_outlier_membership, dtype=np.float32).copy()
    lesion_relation = np.asarray(initial_lesion_relation, dtype=np.float32).copy()
    vent = np.asarray(ventricle_distances, dtype=np.int16)
    gm = np.asarray(gray_matter_distances, dtype=np.int16)
    brainstem = np.asarray(brainstem_distances, dtype=np.int16)
    if current.shape != (*spatial, 13) or shapes.shape != current.shape or direct.shape != current.shape:
        raise ValueError("class membership, shape, and relation state must match the 13-class spatial grid")
    if any(array.shape != spatial for array in (lesion, black_hole, outlier, lesion_relation, vent, gm, brainstem)):
        raise ValueError("auxiliary membership, relation, and distance maps must match the spatial grid")
    if not np.isfinite(data).all() or not np.isfinite(shapes).all() or not np.isfinite(current).all():
        raise ValueError("images, transformed shapes, and initial memberships must be finite")
    for name, array in (("lesion", lesion), ("black_hole", black_hole),
                        ("outlier", outlier), ("lesion_relation", lesion_relation)):
        if not np.isfinite(array).all():
            raise ValueError(f"{name} map must be finite")

    class_energy = np.zeros_like(current)
    lesion_energy = np.zeros(spatial, dtype=np.float32)
    black_hole_energy = np.zeros(spatial, dtype=np.float32)
    denominator = np.zeros(spatial, dtype=np.float32)
    nx, ny, nz = spatial
    # LesionToadSegmentation computes centroid similarities once before its
    # ordered voxel sweep. Keeping this invariant source state outside the
    # sweep avoids rebuilding the same 13x13 matrices at every voxel.
    similarities = compute_lesion_toads_similarity(
        class_centroids,
        class_stddev,
        lesion_centroid,
        lesion_stddev,
        channel_weights,
        tuple(kind.lower() for kind in STANDARD13_TYPES),
        STANDARD13_NAMES,
        prior_scale=float(STANDARD13_PRIOR_SCALE),
        t1_channel=0,
        black_hole_similarity=np.zeros(13, dtype=np.float32),
        finite_zero_variance=not source_zero_variance_nan,
    )

    # Java uses x outermost, then y, then z, mutating maps before the next
    # coordinate reads its six neighbors.
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                faces = np.stack((current[x + 1, y, z], current[x - 1, y, z],
                                  current[x, y + 1, z], current[x, y - 1, z],
                                  current[x, y, z + 1], current[x, y, z - 1])).astype(np.float32)
                lesion_faces = np.asarray((lesion[x + 1, y, z], lesion[x - 1, y, z],
                                           lesion[x, y + 1, z], lesion[x, y - 1, z],
                                           lesion[x, y, z + 1], lesion[x, y, z - 1]), dtype=np.float32)
                black_hole_faces = np.asarray((black_hole[x + 1, y, z], black_hole[x - 1, y, z],
                                               black_hole[x, y + 1, z], black_hole[x, y - 1, z],
                                               black_hole[x, y, z + 1], black_hole[x, y, z - 1]), dtype=np.float32)
                voxel = compute_lesion_toads_standard13_voxel(
                    data[x, y, z], class_centroids, class_stddev, lesion_centroid, lesion_stddev,
                    channel_weights, shapes[x, y, z], shape_priors, smoothing_priors,
                    direct[x, y, z], current[x, y, z], lesion[x, y, z], black_hole[x, y, z],
                    outlier[x, y, z], lesion_relation[x, y, z], faces, lesion_faces,
                    black_hole_faces, outlier_energy=outlier_energy,
                    lesion_power_exponent=lesion_power_exponent,
                    ventricle_distance=int(vent[x, y, z]),
                    gray_matter_distance=int(gm[x, y, z]),
                    brainstem_distance=int(brainstem[x, y, z]), spread=spread,
                    gm_flag=gm_flag, source_zero_variance_nan=source_zero_variance_nan,
                    intensity_field=None if fields is None else fields[x, y, z],
                    precomputed_similarities=similarities,
                )
                current[x, y, z] = voxel.class_memberships
                lesion[x, y, z] = voxel.lesion_membership
                black_hole[x, y, z] = voxel.black_hole_membership
                outlier[x, y, z] = voxel.outlier_membership
                lesion_relation[x, y, z] = voxel.lesion_relation
                class_energy[x, y, z] = voxel.class_energies
                lesion_energy[x, y, z] = voxel.lesion_energy
                black_hole_energy[x, y, z] = voxel.black_hole_energy
                denominator[x, y, z] = voxel.denominator

    objective_result = reduce_lesion_toads_objective(
        current, class_energy, denominator, lesion, black_hole,
        lesion_energy, black_hole_energy,
    )
    return Standard13VolumeResult(
        class_memberships=current,
        lesion_membership=lesion,
        black_hole_membership=black_hole,
        outlier_membership=outlier,
        lesion_relation=lesion_relation,
        denominator=denominator,
        class_energies=class_energy,
        lesion_energy=lesion_energy,
        black_hole_energy=black_hole_energy,
        objective=objective_result.total,
        running_objective=objective_result.running_total,
    )
