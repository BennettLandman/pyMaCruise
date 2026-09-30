"""Selected LesionTOADS class data-energy accumulators.

Paper context: Shiee et al. (2010), ``docs/papers/lesion-toads.md``, describes
statistical intensity terms for lesion-aware classification. P39v validates
the selected field-multiplied image branch alongside the unfielded terms.
These accumulators use fixed inputs and do not establish the complete paper
energy model or iteration.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.toads_math import toads_powertable_lookup


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
    flair_channel: int = -1,
    use_lesion_weight: bool = False,
    lesion_power_exponent: float | None = None,
    cached_lesion_power_exponent: float | None = None,
    structure_name: str | None = None,
    t1_channel: int = -1,
    has_mprage: bool = False,
    csf_power_exponent: float | None = None,
    cached_csf_power_exponent: float | None = None,
) -> LesionToadsDataEnergies:
    """Compute one selected unfielded class's weighted squared residuals.

    The implementation matches the selected Java `OBJECT`, `LESION`, and
    `BLACKHOLE` loops before power weighting, inversion, smoothing, atlas
    priors, and membership normalization. With ``use_field=True``, each image
    channel is multiplied by its field value before centroid residuals, as in
    the Java branch. `obj` is accepted as the repository descriptor's name
    for the source `OBJECT` type. Inputs are one voxel's channel vectors,
    already normalized/scaled by the caller. For exact unfielded
    `Sulcal-CSF`, the selected Java source adds its T1 special term and then
    the generic weighted multichannel sum; with an intensity field, the
    special T1 term is exclusive because that branch has an active `else`.
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
    if not isinstance(use_lesion_weight, (bool, np.bool_)):
        raise ValueError("use_lesion_weight must be boolean")
    if not isinstance(flair_channel, (int, np.integer)) or not -1 <= flair_channel < image.size:
        raise ValueError("flair_channel must be -1 or a valid channel index")
    weighted_lesion_energy = bool(
        object_type == "lesion" and use_lesion_weight and flair_channel != -1
    )
    if weighted_lesion_energy:
        if lesion_power_exponent is None or not np.isfinite(lesion_power_exponent):
            raise ValueError("a finite lesion_power_exponent is required for FLAIR lesion weighting")
        cached_exponent = (
            lesion_power_exponent
            if cached_lesion_power_exponent is None
            else cached_lesion_power_exponent
        )
        if not np.isfinite(cached_exponent):
            raise ValueError("cached_lesion_power_exponent must be finite")
    else:
        cached_exponent = None
    if not isinstance(t1_channel, (int, np.integer)) or not -1 <= t1_channel < image.size:
        raise ValueError("t1_channel must be -1 or a valid channel index")
    if not isinstance(has_mprage, (bool, np.bool_)):
        raise ValueError("has_mprage must be boolean")
    sulcal_csf_t1 = bool(
        object_type in {"obj", "object"}
        and structure_name == "Sulcal-CSF"
        and t1_channel != -1
    )
    if sulcal_csf_t1 and has_mprage:
        if csf_power_exponent is None or not np.isfinite(csf_power_exponent):
            raise ValueError("a finite csf_power_exponent is required for Sulcal-CSF MPRAGE weighting")
        cached_csf_exponent = (
            csf_power_exponent
            if cached_csf_power_exponent is None
            else cached_csf_power_exponent
        )
        if not np.isfinite(cached_csf_exponent):
            raise ValueError("cached_csf_power_exponent must be finite")
    else:
        cached_csf_exponent = None

    class_energy = np.float32(0.0)
    lesion_energy = np.float32(0.0)
    black_hole_energy = np.float32(0.0)
    # In the selected unfielded source, the closing `else` before the
    # multichannel loop is inside a block comment. Sulcal-CSF therefore adds
    # its special T1 term first and then falls through to the ordinary
    # weighted all-channel sum. In the fielded branch the `else` is active,
    # so the Sulcal-CSF T1 term replaces that generic sum.
    if sulcal_csf_t1 and not use_field:
        difference = np.float32(values[t1_channel] - centroid[t1_channel])
        squared_difference = np.float32(difference * difference)
        if has_mprage:
            lookup = toads_powertable_lookup(
                float(squared_difference),
                float(csf_power_exponent),
                cached_exponent=float(cached_csf_exponent),
                minimum=0.0,
                maximum=4.0,
                step=float(np.float32(1.0e-6)),
            ).item()
            # Java's float += double lookup promotes to double and narrows
            # after the special term, before the following float loop.
            class_energy = np.float32(float(class_energy) + lookup)
        else:
            class_energy = squared_difference

    for channel in range(image.size):
        if not (sulcal_csf_t1 and use_field):
            difference = np.float32(values[channel] - centroid[channel])
            weighted_difference = np.float32(weights[channel] * difference)
            contribution = np.float32(weighted_difference * difference)
            class_energy = np.float32(class_energy + contribution)
        if object_type == "lesion" and not weighted_lesion_energy:
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

    if sulcal_csf_t1 and use_field:
        difference = np.float32(values[t1_channel] - centroid[t1_channel])
        squared_difference = np.float32(difference * difference)
        if has_mprage:
            lookup = toads_powertable_lookup(
                float(squared_difference),
                float(csf_power_exponent),
                cached_exponent=float(cached_csf_exponent),
                minimum=0.0,
                maximum=4.0,
                step=float(np.float32(1.0e-6)),
            ).item()
            class_energy = np.float32(float(class_energy) + lookup)
        else:
            class_energy = squared_difference

    if weighted_lesion_energy:
        lesion_difference = np.float32(values[flair_channel] - lesion[flair_channel])
        squared_difference = np.float32(lesion_difference * lesion_difference)
        lookup = toads_powertable_lookup(
            float(squared_difference),
            float(lesion_power_exponent),
            cached_exponent=float(cached_exponent),
            minimum=0.0,
            maximum=4.0,
            step=float(np.float32(1.0e-6)),
        ).item()
        # Java's float += double lookup promotes the accumulator, adds in
        # double precision, then narrows the assigned result back to float.
        lesion_energy = np.float32(float(lesion_energy) + lookup)

    return LesionToadsDataEnergies(class_energy, lesion_energy, black_hole_energy)
