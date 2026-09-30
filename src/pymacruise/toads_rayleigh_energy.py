"""Fixed-input volume data-energy adapter for TOADS RAYLEIGH classes."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_rayleigh import toads_rayleigh_distance


def toads_rayleigh_intensity_energy(
    image_values: ArrayLike,
    centroids: ArrayLike,
    modality_weights: ArrayLike,
    class_variance_weights: ArrayLike,
    image_variances: ArrayLike,
    *,
    intensity_field: ArrayLike | None = None,
    psf_bias: ArrayLike | None = None,
) -> NDArray[np.float32]:
    """Compute the selected interior RAYLEIGH data-energy volume.

    For each interior voxel, class, then modality, this applies the selected
    source's field correction and PSF subtraction, calls
    :func:`toads_rayleigh_distance`, multiplies modality × class-variance ×
    image-variance × distance in Java float32 order, and accumulates into the
    class energy. Faces remain zero because the source exact-membership sweep
    visits only interior voxels. All statistical inputs are fixed inputs;
    this function estimates no centroids or variances and performs no
    membership update.

    Image-like arrays use ``(x,y,z,modality)`` order and centroids/variance
    weights use ``(modality,class)`` order. The caller supplies RAS-aligned
    grids, matching the package's array-kernel convention.
    """

    images = np.asarray(image_values, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    modality = np.asarray(modality_weights, dtype=np.float32)
    class_weights = np.asarray(class_variance_weights, dtype=np.float32)
    image_var = np.asarray(image_variances, dtype=np.float32)
    if images.ndim != 4 or images.shape[-1] < 1:
        raise ValueError("image_values must have shape (x,y,z,modality)")
    channels = images.shape[-1]
    if centers.ndim != 2 or centers.shape[0] != channels or centers.shape[1] < 1:
        raise ValueError("centroids must have shape (modality,class)")
    classes = centers.shape[1]
    if modality.shape != (channels,) or image_var.shape != (channels,):
        raise ValueError("modality_weights and image_variances must match modality count")
    if class_weights.shape != (channels, classes):
        raise ValueError("class_variance_weights must match centroid dimensions")

    if intensity_field is None:
        corrected = images
    else:
        field = np.asarray(intensity_field, dtype=np.float32)
        if field.shape != images.shape:
            raise ValueError("intensity_field must match image_values")
        corrected = np.multiply(field, images, dtype=np.float32)
    if psf_bias is None:
        bias = np.zeros_like(images, dtype=np.float32)
    else:
        bias = np.asarray(psf_bias, dtype=np.float32)
        if bias.shape != images.shape:
            raise ValueError("psf_bias must match image_values")
    corrected = np.subtract(corrected, bias, dtype=np.float32)

    energy = np.zeros(images.shape[:3] + (classes,), dtype=np.float32)
    if any(size < 3 for size in images.shape[:3]):
        return energy
    interior = (slice(1, -1), slice(1, -1), slice(1, -1))
    for class_index in range(classes):
        total = np.zeros(images.shape[:3], dtype=np.float32)
        for channel in range(channels):
            distance = toads_rayleigh_distance(
                corrected[interior + (channel,)], centers[channel, class_index]
            )
            weight = np.float32(modality[channel] * class_weights[channel, class_index])
            weight = np.float32(weight * image_var[channel])
            contribution = np.multiply(weight, distance, dtype=np.float32)
            total[interior] = np.add(total[interior], contribution, dtype=np.float32)
        energy[..., class_index] = total
    return energy
