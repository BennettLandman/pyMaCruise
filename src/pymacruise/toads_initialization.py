"""Source-faithful TOADS intensity and exact-membership initialization.

The formulas in this module follow the selected ``MultipleToadPreprocess``
implementation. Arrays use ``(channel, x, y, z)`` order, matching the Java
``float[channel][x][y][z]`` representation and the selected wrapper's exact
membership path. Values are evaluated as float32 in source iteration order.

Source anchor: ``toads-cruise/src/edu/jhmi/rad/medic/methods/MultipleToadPreprocess.java``
(March 2007; Pierre-Louis Bazin). The selected wrapper/driver path and revision
are recorded in ``docs/toads-port-plan.md`` and ``docs/assets-and-dependencies.md``.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


_ZERO = np.float32(1.0e-20)
_INF = np.float32(1.0e20)


def initial_intensity_centroids(
    low: ArrayLike,
    high: ArrayLike,
    priors: ArrayLike,
) -> NDArray[np.float32]:
    """Map atlas intensity priors into the robust-normalized image range.

    This ports ``MultipleToadPreprocess.initialCentroids``. Since image
    normalization divides by ``high-low`` without subtracting ``low``, each
    class centroid is ``low/(high-low) + prior``. ``low`` and ``high`` are
    one value per input channel; ``priors`` has shape ``(channels, classes)``.
    """
    low_values = np.asarray(low, dtype=np.float32)
    high_values = np.asarray(high, dtype=np.float32)
    prior_values = np.asarray(priors, dtype=np.float32)
    if low_values.ndim != 1 or high_values.shape != low_values.shape:
        raise ValueError("low and high must be matching one-dimensional channel arrays")
    if prior_values.ndim != 2 or prior_values.shape[0] != low_values.size:
        raise ValueError("priors must have shape (channels, classes)")
    if prior_values.shape[1] < 1:
        raise ValueError("priors must contain at least one class")
    if not (
        np.isfinite(low_values).all()
        and np.isfinite(high_values).all()
        and np.isfinite(prior_values).all()
    ):
        raise ValueError("intensity bounds and priors must be finite")
    span = np.subtract(high_values, low_values, dtype=np.float32)
    if np.any(span <= 0):
        raise ValueError("every channel must have a positive robust intensity span")
    offset = np.divide(low_values, span, dtype=np.float32)
    return np.add(offset[:, None], prior_values, dtype=np.float32)


def initial_exact_memberships(
    images: ArrayLike,
    centroids: ArrayLike,
    *,
    crop_border: int = 10,
) -> NDArray[np.float32]:
    """Compute exact reciprocal-squared-distance initial memberships.

    The selected ``AlgorithmMultipleToads`` wrapper sets ``useApprox=false``.
    ``MultipleToadPreprocess.initialExactMemberships`` writes only inside the
    crop: from ``crop_border`` through (but excluding) the last unpadded source
    voxel on each axis. The padded exterior and the crop's last source plane
    remain zero, as in its newly allocated Java output array.

    At each written voxel, squared differences are accumulated in channel
    order. Distances at or below ``1e-20f`` map to ``1e20f``. Each reciprocal
    is multiplied by the number of classes whose first-channel centroid is
    exactly equal, then normalized across classes.
    """
    image_values = np.asarray(images, dtype=np.float32)
    centroid_values = np.asarray(centroids, dtype=np.float32)
    if image_values.ndim != 4 or any(size == 0 for size in image_values.shape):
        raise ValueError("images must have shape (channels, x, y, z) with nonempty axes")
    if centroid_values.ndim != 2 or centroid_values.shape[0] != image_values.shape[0]:
        raise ValueError("centroids must have shape (channels, classes)")
    if centroid_values.shape[1] == 0:
        raise ValueError("centroids must contain at least one class")
    if crop_border < 0:
        raise ValueError("crop_border must be non-negative")
    if not np.isfinite(image_values).all() or not np.isfinite(centroid_values).all():
        raise ValueError("images and centroids must be finite")

    bounds = tuple(
        slice(crop_border, size - crop_border - 1)
        for size in image_values.shape[1:]
    )
    if any(item.stop <= item.start for item in bounds):
        raise ValueError("image is too small for the requested TOADS crop border")

    classes = centroid_values.shape[1]
    memberships = np.zeros(
        image_values.shape[1:] + (classes,), dtype=np.float32
    )
    multiplicity = np.zeros(classes, dtype=np.float32)
    for k in range(classes):
        for m in range(classes):
            if centroid_values[0, m] == centroid_values[0, k]:
                multiplicity[k] = np.float32(multiplicity[k] + np.float32(1.0))

    for x in range(bounds[0].start, bounds[0].stop):
        for y in range(bounds[1].start, bounds[1].stop):
            for z in range(bounds[2].start, bounds[2].stop):
                numerators = np.empty(classes, dtype=np.float32)
                denominator = np.float32(0.0)
                for k in range(classes):
                    distance = np.float32(0.0)
                    for channel in range(image_values.shape[0]):
                        difference = np.float32(
                            image_values[channel, x, y, z]
                            - centroid_values[channel, k]
                        )
                        squared = np.float32(difference * difference)
                        distance = np.float32(distance + squared)
                    if distance > _ZERO:
                        numerator = np.float32(np.float32(1.0) / distance)
                    else:
                        numerator = _INF
                    numerators[k] = numerator
                    denominator = np.float32(denominator + numerator)

                if denominator > np.float32(0.0):
                    for k in range(classes):
                        weighted = np.float32(multiplicity[k] * numerators[k])
                        memberships[x, y, z, k] = np.float32(weighted / denominator)

    return memberships
