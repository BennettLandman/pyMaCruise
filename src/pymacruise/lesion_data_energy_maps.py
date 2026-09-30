"""Volume-level raw data-energy maps for selected LesionTOADS classes.

This adapter composes the already validated one-voxel energy helper over the
selected ``computeMemberships`` interior scan. It stops before smoothing,
shape priors, inversion, membership changes, or objective accumulation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.lesion_data_terms import compute_lesion_toads_unfielded_data_energies


@dataclass(frozen=True)
class LesionToadsDataEnergyMaps:
    """Raw class and shared auxiliary energy maps, all stored as float32."""

    class_energy: NDArray[np.float32]
    lesion_energy: NDArray[np.float32]
    black_hole_energy: NDArray[np.float32]


def compute_lesion_toads_data_energy_maps(
    images: ArrayLike,
    class_centroids: ArrayLike,
    lesion_centroid: ArrayLike,
    black_hole_centroid: ArrayLike,
    channel_weights: ArrayLike,
    class_types: Sequence[str],
) -> LesionToadsDataEnergyMaps:
    """Compute raw unfielded OBJECT/LESION/BLACKHOLE energies over a volume.

    ``images`` is ``(x, y, z, channels)`` and ``class_centroids`` is
    ``(classes, channels)``. Only the selected Java interior (indices
    ``1..n-2`` on every spatial axis) is evaluated; energy maps retain zeros
    on the untouched one-voxel faces. Auxiliary LESION and BLACKHOLE energies
    share one per-voxel slot in class order: each LESION or BLACKHOLE class
    overwrites lesion energy, while each BLACKHOLE class also overwrites the
    black-hole energy slot. A later LESION therefore leaves the most recent
    black-hole value intact.

    This fixed-input adapter does not model topology/relation skips, fields,
    Sulcal-CSF or PowerTable behavior, FLAIR weighting, spatial energy terms,
    inversion, membership normalization, or objective reduction.
    """
    values = np.asarray(images, dtype=np.float32)
    centroids = np.asarray(class_centroids, dtype=np.float32)
    lesion = np.asarray(lesion_centroid, dtype=np.float32)
    black_hole = np.asarray(black_hole_centroid, dtype=np.float32)
    weights = np.asarray(channel_weights, dtype=np.float32)

    if values.ndim != 4 or min(values.shape[:3], default=0) < 3 or values.shape[3] < 1:
        raise ValueError("images must have shape (x>=3, y>=3, z>=3, channels>=1)")
    if centroids.ndim != 2 or centroids.shape[0] < 1:
        raise ValueError("class_centroids must have shape (classes>=1, channels)")
    channels = values.shape[3]
    if centroids.shape[1] != channels:
        raise ValueError("class_centroids channel count must match images")
    for name, vector in (
        ("lesion_centroid", lesion),
        ("black_hole_centroid", black_hole),
        ("channel_weights", weights),
    ):
        if vector.shape != (channels,):
            raise ValueError(f"{name} must have shape (channels,)")
    if len(class_types) != centroids.shape[0]:
        raise ValueError("class_types length must match class_centroids")
    if not np.isfinite(values).all():
        raise ValueError("images must contain only finite values")
    for name, array in (
        ("class_centroids", centroids),
        ("lesion_centroid", lesion),
        ("black_hole_centroid", black_hole),
        ("channel_weights", weights),
    ):
        if not np.isfinite(array).all():
            raise ValueError(f"{name} must contain only finite values")

    normalized_types: list[str] = []
    aliases = {"obj": "object", "object": "object", "lesion": "lesion", "blackhole": "blackhole"}
    for item in class_types:
        if not isinstance(item, str) or item.lower() not in aliases:
            raise ValueError("class_types may contain only OBJECT, LESION, or BLACKHOLE")
        normalized_types.append(aliases[item.lower()])

    shape = values.shape[:3]
    class_energy = np.zeros((*shape, len(normalized_types)), dtype=np.float32)
    lesion_energy = np.zeros(shape, dtype=np.float32)
    black_hole_energy = np.zeros(shape, dtype=np.float32)

    # The Java source scans x, then y, then z and evaluates classes in their
    # declared order. Keep this explicit order because auxiliary energy slots
    # are shared and overwritten by subsequent lesion-related classes.
    for x in range(1, shape[0] - 1):
        for y in range(1, shape[1] - 1):
            for z in range(1, shape[2] - 1):
                for k, object_type in enumerate(normalized_types):
                    result = compute_lesion_toads_unfielded_data_energies(
                        values[x, y, z, :],
                        centroids[k],
                        lesion,
                        black_hole,
                        weights,
                        object_type=object_type,
                    )
                    class_energy[x, y, z, k] = result.class_energy
                    if object_type in {"lesion", "blackhole"}:
                        lesion_energy[x, y, z] = result.lesion_energy
                    if object_type == "blackhole":
                        black_hole_energy[x, y, z] = result.black_hole_energy

    return LesionToadsDataEnergyMaps(class_energy, lesion_energy, black_hole_energy)
