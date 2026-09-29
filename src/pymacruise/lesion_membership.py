"""Selected scalar membership conversions from LesionTOADS updates."""

from __future__ import annotations

import numpy as np


def convert_lesion_toads_class_membership(
    data_energy: float,
    shape_prior: float,
    *,
    relation_type: str,
    zero: float = 1.0e-10,
    infinity: float = 1.0e10,
) -> np.float32:
    """Invert one class energy and apply its direct/prior relation rule.

    This reproduces the selected `computeMemberships` class branch before
    denominator accumulation and voxelwise membership normalization. Values
    use Java float32 semantics. Energies at or below the strict ``ZERO``
    threshold map to the selected finite ``INF`` sentinel.
    """
    energy = np.float32(data_energy)
    shape = np.float32(shape_prior)
    zero_value = np.float32(zero)
    infinity_value = np.float32(infinity)
    if not np.isfinite(energy):
        raise ValueError("data_energy must be finite")
    if not np.isfinite(shape):
        raise ValueError("shape_prior must be finite")
    if not np.isfinite(zero_value) or zero_value < np.float32(0.0):
        raise ValueError("zero must be finite and nonnegative")
    if not np.isfinite(infinity_value) or infinity_value < np.float32(0.0):
        raise ValueError("infinity sentinel must be finite and nonnegative")
    if relation_type not in {"direct", "prior"}:
        raise ValueError("relation_type must be 'direct' or 'prior'")

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        membership = (
            np.float32(np.float32(1.0) / energy)
            if energy > zero_value
            else infinity_value
        )
    if relation_type == "prior":
        doubled_shape = np.float32(np.float32(2.0) * shape)
        # VandyCRUISE Numerics.min(a,b) returns b on ties and when b is NaN.
        shape_factor = np.float32(1.0) if np.float32(1.0) < doubled_shape else doubled_shape
        membership = np.float32(shape_factor * membership)
    return membership
