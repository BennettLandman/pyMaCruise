"""Selected unfielded LesionTOADS class data-energy accumulators."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class LesionToadsDataEnergies:
    """Raw class, lesion, and black-hole data terms before inversion."""

    class_energy: np.float32
    lesion_energy: np.float32
    black_hole_energy: np.float32


def compute_lesion_toads_unfielded_data_energies(
    image_values: ArrayLike,
    class_centroid: ArrayLike,
    lesion_centroid: ArrayLike,
    black_hole_centroid: ArrayLike,
    channel_weight: ArrayLike,
    *,
    object_type: str,
    use_field: bool = False,
    intensity_field: ArrayLike | None = None,
) -> LesionToadsDataEnergies:
    """Compute one selected unfielded class's weighted squared residuals.

    The implementation matches the selected Java `OBJECT`, `LESION`, and
    `BLACKHOLE` loops before power weighting, inversion, smoothing, atlas
    priors, and membership normalization. With ``use_field=True``, each image
    channel is multiplied by its field value before centroid residuals, as in
    the Java branch. `obj` is accepted as the repository descriptor's name
    for the source `OBJECT` type. Inputs are one voxel's channel vectors,
    already normalized/scaled by the caller.
    """
    image = np.asarray(image_values, dtype=np.float32)
    centroid = np.asarray(class_centroid, dtype=np.float32)
    lesion = np.asarray(lesion_centroid, dtype=np.float32)
    black_hole = np.asarray(black_hole_centroid, dtype=np.float32)
    weights = np.asarray(channel_weight, dtype=np.float32)
    if not isinstance(use_field, (bool, np.bool_)):
        raise ValueError("use_field must be boolean")
    if use_field:
        if intensity_field is None:
            raise ValueError("intensity_field is required when use_field is true")
        field = np.asarray(intensity_field, dtype=np.float32)
        if field.shape != image.shape:
            raise ValueError("intensity_field must match image channel count")
        if not np.isfinite(field).all():
            raise ValueError("intensity_field must contain only finite values")
        values = np.asarray(field * image, dtype=np.float32)
    else:
        values = image
    if image.ndim != 1 or image.size == 0:
        raise ValueError("image_values must be a nonempty channel vector")
    for name, array in (
        ("class_centroid", centroid),
        ("lesion_centroid", lesion),
        ("black_hole_centroid", black_hole),
        ("channel_weight", weights),
    ):
        if array.shape != image.shape:
            raise ValueError(f"{name} must match image channel count")
        if not np.isfinite(array).all():
            raise ValueError(f"{name} must contain only finite values")
    if not np.isfinite(image).all():
        raise ValueError("image_values must contain only finite values")
    if object_type not in {"obj", "object", "lesion", "blackhole"}:
        raise ValueError("object_type must be obj/object, lesion, or blackhole")

    class_energy = np.float32(0.0)
    lesion_energy = np.float32(0.0)
    black_hole_energy = np.float32(0.0)
    for channel in range(image.size):
        difference = np.float32(values[channel] - centroid[channel])
        weighted_difference = np.float32(weights[channel] * difference)
        contribution = np.float32(weighted_difference * difference)
        class_energy = np.float32(class_energy + contribution)
        if object_type == "lesion":
            lesion_difference = np.float32(values[channel] - lesion[channel])
            lesion_weighted_difference = np.float32(weights[channel] * lesion_difference)
            lesion_contribution = np.float32(lesion_weighted_difference * lesion_difference)
            lesion_energy = np.float32(lesion_energy + lesion_contribution)
        elif object_type == "blackhole":
            lesion_difference = np.float32(values[channel] - lesion[channel])
            lesion_weighted_difference = np.float32(weights[channel] * lesion_difference)
            lesion_contribution = np.float32(lesion_weighted_difference * lesion_difference)
            lesion_energy = np.float32(lesion_energy + lesion_contribution)
            black_hole_difference = np.float32(values[channel] - black_hole[channel])
            black_hole_weighted_difference = np.float32(weights[channel] * black_hole_difference)
            black_hole_contribution = np.float32(black_hole_weighted_difference * black_hole_difference)
            black_hole_energy = np.float32(black_hole_energy + black_hole_contribution)

    return LesionToadsDataEnergies(class_energy, lesion_energy, black_hole_energy)
