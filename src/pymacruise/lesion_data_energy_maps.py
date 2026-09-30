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
    *,
    use_field: bool = False,
    intensity_fields: ArrayLike | None = None,
    structure_names: Sequence[str | None] | None = None,
    t1_channel: int = -1,
    has_mprage: bool = False,
    csf_power_exponent: float | None = None,
    cached_csf_power_exponent: float | None = None,
    flair_channel: int = -1,
    use_lesion_weight: bool = False,
    lesion_power_exponent: float | None = None,
    cached_lesion_power_exponent: float | None = None,
) -> LesionToadsDataEnergyMaps:
    """Compute raw OBJECT/LESION/BLACKHOLE energies over a volume.

    ``images`` is ``(x, y, z, channels)`` and ``class_centroids`` is
    ``(classes, channels)``. Only the selected Java interior (indices
    ``1..n-2`` on every spatial axis) is evaluated; energy maps retain zeros
    on the untouched one-voxel faces. Auxiliary LESION and BLACKHOLE energies
    share one per-voxel slot in class order: each LESION or BLACKHOLE class
    overwrites lesion energy, while each BLACKHOLE class also overwrites the
    black-hole energy slot. A later LESION therefore leaves the most recent
    black-hole value intact.

    Per-channel ``intensity_fields`` are applied before residuals when
    ``use_field`` is true. Exact-name OBJECT ``Sulcal-CSF`` routes through the
    selected T1/MPRAGE branch in the scalar helper; when unfielded, that helper
    also preserves the source's subsequent generic multichannel accumulation.
    LESION's optional FLAIR PowerTable term is likewise delegated to the
    validated scalar helper. This adapter still excludes topology/relation
    skips, spatial energy terms, inversion, membership normalization, and
    objective reduction.
    """
    values = np.asarray(images, dtype=np.float32)
    centroids = np.asarray(class_centroids, dtype=np.float32)
    lesion = np.asarray(lesion_centroid, dtype=np.float32)
    black_hole = np.asarray(black_hole_centroid, dtype=np.float32)
    weights = np.asarray(channel_weights, dtype=np.float32)

    if not isinstance(use_field, (bool, np.bool_)):
        raise ValueError("use_field must be boolean")
    if not isinstance(has_mprage, (bool, np.bool_)):
        raise ValueError("has_mprage must be boolean")
    if not isinstance(use_lesion_weight, (bool, np.bool_)):
        raise ValueError("use_lesion_weight must be boolean")

    if values.ndim != 4 or min(values.shape[:3], default=0) < 3 or values.shape[3] < 1:
        raise ValueError("images must have shape (x>=3, y>=3, z>=3, channels>=1)")
    if centroids.ndim != 2 or centroids.shape[0] < 1:
        raise ValueError("class_centroids must have shape (classes>=1, channels)")
    channels = values.shape[3]
    if centroids.shape[1] != channels:
        raise ValueError("class_centroids channel count must match images")
    if intensity_fields is None:
        fields = None
    else:
        fields = np.asarray(intensity_fields, dtype=np.float32)
        if fields.shape != values.shape:
            raise ValueError("intensity_fields must match images shape")
        if not np.isfinite(fields).all():
            raise ValueError("intensity_fields must contain only finite values")
    if use_field and fields is None:
        raise ValueError("intensity_fields are required when use_field is true")
    if structure_names is None:
        names: tuple[str | None, ...] = (None,) * centroids.shape[0]
    else:
        names = tuple(structure_names)
        if len(names) != centroids.shape[0] or any(name is not None and not isinstance(name, str) for name in names):
            raise ValueError("structure_names must contain one string or None per class")
    for name, channel in (("t1_channel", t1_channel), ("flair_channel", flair_channel)):
        if isinstance(channel, (bool, np.bool_)) or not isinstance(channel, (int, np.integer)) or not -1 <= int(channel) < channels:
            raise ValueError(f"{name} must be -1 or a valid channel index")
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
                        use_field=use_field,
                        intensity_field=None if fields is None else fields[x, y, z, :],
                        flair_channel=flair_channel,
                        use_lesion_weight=use_lesion_weight,
                        lesion_power_exponent=lesion_power_exponent,
                        cached_lesion_power_exponent=cached_lesion_power_exponent,
                        structure_name=names[k],
                        t1_channel=t1_channel,
                        has_mprage=has_mprage,
                        csf_power_exponent=csf_power_exponent,
                        cached_csf_power_exponent=cached_csf_power_exponent,
                    )
                    class_energy[x, y, z, k] = result.class_energy
                    if object_type in {"lesion", "blackhole"}:
                        lesion_energy[x, y, z] = result.lesion_energy
                    if object_type == "blackhole":
                        black_hole_energy[x, y, z] = result.black_hole_energy

    return LesionToadsDataEnergyMaps(class_energy, lesion_energy, black_hole_energy)
