"""Data terms for every active class family in the bundled LesionTOADS atlas.

The selected VandyCRUISE constructor maps the 13-class descriptor to OBJECT,
MASK, OPTIMIZED_CSF, and LESION. In ``computeMemberships`` the optimized-CSF
formula is commented out, so those classes reach the source's fallback scalar
outlier term. This module makes that active dispatch explicit. It only builds
raw data terms; smoothing, priors, memberships, and convergence are separate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_data_terms import compute_lesion_toads_unfielded_data_energies


@dataclass(frozen=True)
class LesionToadsActiveClassEnergies:
    """Raw class, lesion, and black-hole energy maps in descriptor order."""

    class_energy: NDArray[np.float32]
    lesion_energy: NDArray[np.float32]
    black_hole_energy: NDArray[np.float32]


_OUTLIER_FALLBACK_TYPES = {
    "outlier", "rayleigh", "optimized_csf", "optimized_gm", "optimized_wm",
    "optimized_lesion", "optimized_blackhole",
}


def compute_lesion_toads_active_class_energies(
    images: ArrayLike,
    class_centroids: ArrayLike,
    lesion_centroid: ArrayLike,
    black_hole_centroid: ArrayLike,
    channel_weights: ArrayLike,
    class_types: Sequence[str],
    *,
    outlier_distance: float,
    use_field: bool = False,
    intensity_fields: ArrayLike | None = None,
) -> LesionToadsActiveClassEnergies:
    """Build raw data terms for the selected source's active class dispatch.

    Arrays follow source x/y/z/channel order; untouched outer faces stay zero.
    ``MASK`` uses binary mismatch distance to zero. The active source branch
    for optimized types falls through to its scalar ``outlier`` assignment
    because the optimized formulas in this method are commented out.
    """
    values = np.asarray(images, dtype=np.float32)
    centroids = np.asarray(class_centroids, dtype=np.float32)
    lesion = np.asarray(lesion_centroid, dtype=np.float32)
    black_hole = np.asarray(black_hole_centroid, dtype=np.float32)
    weights = np.asarray(channel_weights, dtype=np.float32)
    if values.ndim != 4 or min(values.shape[:3], default=0) < 3 or values.shape[3] < 1:
        raise ValueError("images must have shape (x>=3,y>=3,z>=3,channels>=1)")
    if centroids.shape != (len(class_types), values.shape[3]):
        raise ValueError("class_centroids must match class_types and image channels")
    if lesion.shape != (values.shape[3],) or black_hole.shape != lesion.shape or weights.shape != lesion.shape:
        raise ValueError("centroids and channel_weights must match image channels")
    normalized = tuple(str(kind).lower() for kind in class_types)
    supported = {"obj", "object", "lesion", "blackhole", "mask"} | _OUTLIER_FALLBACK_TYPES
    if not normalized or any(kind not in supported for kind in normalized):
        raise ValueError("class_types contain a type not handled by active VandyCRUISE dispatch")
    if not isinstance(use_field, (bool, np.bool_)):
        raise ValueError("use_field must be boolean")
    if not np.isfinite(values).all() or not np.isfinite(centroids).all() or not np.isfinite(lesion).all() or not np.isfinite(black_hole).all() or not np.isfinite(weights).all():
        raise ValueError("image, centroid, and weight inputs must be finite")
    if not np.isfinite(outlier_distance):
        raise ValueError("outlier_distance must be finite")
    fields = None if intensity_fields is None else np.asarray(intensity_fields, dtype=np.float32)
    if fields is not None and (fields.shape != values.shape or not np.isfinite(fields).all()):
        raise ValueError("intensity_fields must be finite and match images")
    if use_field and fields is None:
        raise ValueError("intensity_fields are required when use_field is true")

    shape = values.shape[:3]
    class_energy = np.zeros((*shape, len(normalized)), dtype=np.float32)
    lesion_energy = np.zeros(shape, dtype=np.float32)
    black_energy = np.zeros(shape, dtype=np.float32)
    for x in range(1, shape[0] - 1):
        for y in range(1, shape[1] - 1):
            for z in range(1, shape[2] - 1):
                for k, kind in enumerate(normalized):
                    if kind in _OUTLIER_FALLBACK_TYPES:
                        energy = np.float32(outlier_distance)
                    elif kind == "mask":
                        energy = np.float32(0.0)
                        for channel in range(values.shape[3]):
                            mismatch = np.float32(0.0 if values[x, y, z, channel] == np.float32(0.0) else 1.0)
                            energy = np.float32(energy + np.float32(weights[channel] * mismatch))
                    else:
                        result = compute_lesion_toads_unfielded_data_energies(
                            values[x, y, z], centroids[k], lesion, black_hole, weights,
                            object_type=kind,
                            use_field=use_field,
                            intensity_field=None if fields is None else fields[x, y, z],
                        )
                        energy = result.class_energy
                        if kind == "lesion":
                            lesion_energy[x, y, z] = result.lesion_energy
                        elif kind == "blackhole":
                            lesion_energy[x, y, z] = result.lesion_energy
                            black_energy[x, y, z] = result.black_hole_energy
                    class_energy[x, y, z, k] = energy
    return LesionToadsActiveClassEnergies(class_energy, lesion_energy, black_energy)
