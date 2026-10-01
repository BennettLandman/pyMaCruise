"""LesionTOADS initial centroid handoff into segmentation state.

Ports the active ``LesionToadSegmentation.initCentroids(float[][])`` loop.
The segmentation constructor has already allocated zero-filled current and
prior centroid matrices; the wrapper then copies each channel/class value
from ``LesionToadPreprocess.initialCentroids`` into both matrices.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass
class LesionToadsCentroidHandoffState:
    """The two Java fields written by ``initCentroids``."""

    centroid: NDArray[np.float32]
    prior_centroid: NDArray[np.float32]

    @classmethod
    def empty(cls, channels: int, classes: int) -> "LesionToadsCentroidHandoffState":
        """Allocate the zero-filled arrays created by the Java constructor."""
        if channels < 1 or classes < 1:
            raise ValueError("channels and classes must be positive")
        shape = (int(channels), int(classes))
        return cls(
            centroid=np.zeros(shape, dtype=np.float32),
            prior_centroid=np.zeros(shape, dtype=np.float32),
        )


def initialize_lesion_toads_segmentation_centroids(
    state: LesionToadsCentroidHandoffState,
    initial_centroids: ArrayLike,
) -> None:
    """Copy initial channel/class means to current and prior centroid fields.

    The order and two independent assignments match Java's nested channel,
    class loops. Java float fields are reproduced by storing into float32
    arrays. The source does not retain the input matrix reference.
    """
    values = np.asarray(initial_centroids, dtype=np.float32)
    channels, classes = state.centroid.shape
    for channel in range(channels):
        for class_index in range(classes):
            value = np.float32(values[channel, class_index])
            state.centroid[channel, class_index] = value
            state.prior_centroid[channel, class_index] = value
