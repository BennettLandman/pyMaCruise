"""Fixed-input LesionTOADS FCM membership initialization.

This ports the selected VandyCRUISE ``LesionToadPreprocess.initialMemberships``
loop after preprocessing has produced its padded/cropped normalized image
arrays and tissue centroids. It does not implement image normalization,
cropping, centroid estimation, or atlas registration.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

_ZERO = np.float32(1.0e-20)
_INF = np.float32(1.0e20)


@dataclass(frozen=True)
class LesionInitialMemberships:
    """New membership array and its active padded crop interior."""

    memberships: NDArray[np.float32]
    active_bounds: tuple[tuple[int, int], tuple[int, int], tuple[int, int]]


def lesion_initial_memberships(
    padded_images: ArrayLike,
    centroids: ArrayLike,
    *,
    crop_bounds: tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
    pad_border: int,
) -> LesionInitialMemberships:
    """Compute selected exact FCM memberships in a padded cropped image.

    ``padded_images`` has ``(channels, x, y, z)`` order, follows the caller's
    RAS voxel axes, and is the output of the source crop operation, including
    ``pad_border`` voxels on each face.
    ``crop_bounds`` supplies the source ``(x0i,xNi)``, ``(y0i,yNi)``,
    ``(z0i,zNi)`` pairs; their extents must match the non-padding image region.
    The new result has ``(x,y,z,classes)`` order, is zero outside the active
    crop, and is written at index ``source_coordinate - crop_start + pad``.
    
    At each active voxel, source order is class then channel for squared
    distance and class then class for exact channel-0 centroid multiplicity.
    Distances above ``1e-20f`` are reciprocated; zero/small distances use
    ``1e20f``. Ties in channel-0 centroids multiply a class's reciprocal
    distance by the exact tied-class count before normalization.
    """
    images = np.asarray(padded_images, dtype=np.float32)
    means = np.asarray(centroids, dtype=np.float32)
    if images.ndim != 4 or any(size == 0 for size in images.shape):
        raise ValueError("padded_images must have shape (channels, x, y, z)")
    if means.ndim != 2 or means.shape[0] != images.shape[0] or means.shape[1] == 0:
        raise ValueError("centroids must have shape (channels, nonzero classes)")
    if not np.isfinite(images).all() or not np.isfinite(means).all():
        raise ValueError("padded images and centroids must be finite")
    if isinstance(pad_border, (bool, np.bool_)) or not isinstance(pad_border, (int, np.integer)) or pad_border < 0:
        raise ValueError("pad_border must be a nonnegative integer")
    if len(crop_bounds) != 3:
        raise ValueError("crop_bounds must contain x, y, and z bounds")
    bounds: list[tuple[int, int]] = []
    for axis, pair in enumerate(crop_bounds):
        if len(pair) != 2:
            raise ValueError("each crop bound must be a (start, stop) pair")
        start, stop = pair
        if isinstance(start, (bool, np.bool_)) or not isinstance(start, (int, np.integer)):
            raise ValueError("crop-bound coordinates must be integers")
        if isinstance(stop, (bool, np.bool_)) or not isinstance(stop, (int, np.integer)):
            raise ValueError("crop-bound coordinates must be integers")
        extent = int(stop) - int(start)
        if extent <= 0 or extent + 2 * int(pad_border) != images.shape[axis + 1]:
            raise ValueError("crop extents plus padding must match image dimensions")
        bounds.append((int(start), int(stop)))

    classes = means.shape[1]
    result = np.zeros(images.shape[1:] + (classes,), dtype=np.float32)
    active = tuple((int(pad_border), size - int(pad_border)) for size in images.shape[1:])
    for sx in range(bounds[0][0], bounds[0][1]):
        x = sx - bounds[0][0] + int(pad_border)
        for sy in range(bounds[1][0], bounds[1][1]):
            y = sy - bounds[1][0] + int(pad_border)
            for sz in range(bounds[2][0], bounds[2][1]):
                z = sz - bounds[2][0] + int(pad_border)
                denominator = np.float32(0.0)
                numerators = np.empty(classes, dtype=np.float32)
                multiplicities = np.empty(classes, dtype=np.float32)
                for k in range(classes):
                    distance = np.float32(0.0)
                    for channel in range(images.shape[0]):
                        difference = np.float32(images[channel, x, y, z] - means[channel, k])
                        squared = np.float32(difference * difference)
                        distance = np.float32(distance + squared)
                    numerator = np.float32(np.float32(1.0) / distance) if distance > _ZERO else _INF
                    numerators[k] = numerator
                    denominator = np.float32(denominator + numerator)

                    # The selected method recomputes this inner grouping count
                    # for every class in the same voxel.
                    grouping = np.float32(0.0)
                    for m in range(classes):
                        if means[0, m] == means[0, k]:
                            grouping = np.float32(grouping + np.float32(1.0))
                    multiplicities[k] = grouping

                if denominator > np.float32(0.0):
                    for k in range(classes):
                        weighted = np.float32(multiplicities[k] * numerators[k])
                        result[x, y, z, k] = np.float32(weighted / denominator)
    return LesionInitialMemberships(result, active)  # type: ignore[arg-type]
