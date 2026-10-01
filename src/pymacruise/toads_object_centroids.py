"""Bounded source-faithful ordinary OBJECT centroid update for VandyCRUISE.

This implements only the ordinary ``OBJECT`` numerator/denominator branch
inside ``MultipleToadSegmentation.computeExactCentroids`` and its ``damped``
centroid update. Lesion, outlier, optimized tissue, Rayleigh, linked/prior
modes, and centroid/variance dispatch are outside this adapter.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``, describe
statistical tissue classification using class intensity models. This function
implements only the ordinary-OBJECT centroid update in that broad role; its
summation order and damping details are selected VandyCRUISE semantics, not a
full implementation of the paper's statistical model.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def toads_object_centroid_update(
    images: ArrayLike,
    memberships: ArrayLike,
    previous_centroids: ArrayLike,
    *,
    intensity_field: ArrayLike | None = None,
    segmentation: ArrayLike | None = None,
    psf_ridge: float = 0.0,
    smoothness: float = 0.1,
) -> NDArray[np.float32]:
    """Return one exact-centroid OBJECT accumulation followed by damping.

    Arrays use ``[x,y,z,channel]`` for images/field and ``[x,y,z,class]`` for
    memberships; centroids are ``[channel,class]``. The Java source visits
    only the one-voxel interior, then classes, x, y, z in that nested order.
    A zero denominator in damped mode leaves the incoming centroid unchanged.

    The field, when supplied, multiplies only the center image sample exactly
    as the source does. With ``psf_ridge > 0``, ``segmentation`` is required;
    the original 3x3x3 ``isBoundary`` predicate and boundary Laplacian term
    are applied in their source loop order. The source's separate PSF-spread
    implementation is commented out and is not represented here.
    """
    image = np.asarray(images, dtype=np.float32)
    mem = np.asarray(memberships, dtype=np.float32)
    old = np.asarray(previous_centroids, dtype=np.float32)
    if image.ndim != 4 or mem.ndim != 4 or old.ndim != 2:
        raise ValueError("images, memberships, and previous_centroids must be 4D, 4D, and 2D")
    if image.shape[:3] != mem.shape[:3] or old.shape != (image.shape[3], mem.shape[3]):
        raise ValueError("image, membership, and centroid dimensions must agree")
    if any(size < 3 for size in image.shape[:3]):
        raise ValueError("each spatial dimension needs a one-voxel halo")
    field = None if intensity_field is None else np.asarray(intensity_field, dtype=np.float32)
    if field is not None and field.shape != image.shape:
        raise ValueError("intensity_field must match images")
    ridge = np.float32(psf_ridge)
    if not np.isfinite(ridge) or ridge < np.float32(0.0):
        raise ValueError("psf_ridge must be finite and nonnegative")
    labels = None if segmentation is None else np.asarray(segmentation)
    if ridge > np.float32(0.0):
        if labels is None or labels.shape != image.shape[:3]:
            raise ValueError("segmentation matching the spatial grid is required for PSF ridge")
    elif labels is not None and labels.shape != image.shape[:3]:
        raise ValueError("segmentation must match the image spatial grid")
    blend = np.float32(smoothness)
    if not np.isfinite(blend) or blend < np.float32(0.0) or blend > np.float32(1.0):
        raise ValueError("smoothness must be finite and in [0, 1]")

    classes = mem.shape[3]
    channels = image.shape[3]
    result = old.copy()
    for channel in range(channels):
        numerator = np.zeros(classes, dtype=np.float32)
        denominator = np.zeros(classes, dtype=np.float32)
        for class_index in range(classes):
            for x in range(1, image.shape[0] - 1):
                for y in range(1, image.shape[1] - 1):
                    for z in range(1, image.shape[2] - 1):
                        membership = mem[x, y, z, class_index]
                        sample = (
                            np.float32(field[x, y, z, channel] * image[x, y, z, channel])
                            if field is not None
                            else image[x, y, z, channel]
                        )
                        if ridge > np.float32(0.0) and _is_boundary(labels, x, y, z):
                            sample = np.float32(
                                sample
                                + np.float32(
                                    np.float32(ridge * np.float32(np.float32(28.0) / np.float32(27.0)))
                                    * image[x, y, z, channel]
                                )
                            )
                            for i in range(-1, 2):
                                for j in range(-1, 2):
                                    for k in range(-1, 2):
                                        contribution = np.float32(
                                            np.float32(ridge * np.float32(np.float32(1.0) / np.float32(27.0)))
                                            * image[x + i, y + j, z + k, channel]
                                        )
                                        sample = np.float32(sample - contribution)
                        weighted = np.float32(np.float32(membership * membership) * sample)
                        numerator[class_index] = np.float32(numerator[class_index] + weighted)
                        squared_membership = np.float32(membership * membership)
                        denominator[class_index] = np.float32(
                            denominator[class_index] + squared_membership
                        )
        one_minus_blend = np.float32(np.float32(1.0) - blend)
        for class_index in range(classes):
            if denominator[class_index] > np.float32(0.0):
                undamped = np.float32(
                    np.float32(one_minus_blend * numerator[class_index])
                    / denominator[class_index]
                )
                retained = np.float32(blend * result[channel, class_index])
                result[channel, class_index] = np.float32(undamped + retained)
    return result


def _is_boundary(labels: NDArray, x: int, y: int, z: int) -> bool:
    center = labels[x, y, z]
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                if labels[x + i, y + j, z + k] != center:
                    return True
    return False
