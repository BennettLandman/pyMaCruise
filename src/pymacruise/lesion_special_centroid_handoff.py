"""LesionTOADS lesion and black-hole centroid reference setters.

The active wrapper calls ``setLesionCentroid`` then
``setBlackHoleCentroid`` after the generic centroid handoff. The selected Java
methods assign the supplied ``float[]`` references directly; they do not copy,
convert, or estimate the values. This small state adapter models only those
two assignments.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass
class LesionToadsSpecialCentroidState:
    """The mutable two-vector state changed by the selected Java setters."""

    lesion_centroid: NDArray[np.float32] | None = None
    black_hole_centroid: NDArray[np.float32] | None = None


def set_lesion_centroid_reference(
    state: LesionToadsSpecialCentroidState,
    centroid: NDArray[np.float32],
) -> None:
    """Mirror ``setLesionCentroid(float[])`` by retaining the array reference."""
    state.lesion_centroid = centroid


def set_black_hole_centroid_reference(
    state: LesionToadsSpecialCentroidState,
    centroid: NDArray[np.float32],
) -> None:
    """Mirror ``setBlackHoleCentroid(float[])`` by retaining the array reference."""
    state.black_hole_centroid = centroid


def initialize_lesion_toads_special_centroids(
    state: LesionToadsSpecialCentroidState,
    lesion_centroid: NDArray[np.float32],
    black_hole_centroid: NDArray[np.float32],
) -> None:
    """Apply both setters in the active wrapper's order."""
    set_lesion_centroid_reference(state, lesion_centroid)
    set_black_hole_centroid_reference(state, black_hole_centroid)
