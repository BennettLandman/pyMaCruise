"""Small numerical helpers for source-compatible TOADS likelihoods."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

_SOURCE_STEP = float(np.float32(1e-6))

def toads_powertable_lookup(
    values: ArrayLike,
    exponent: float,
    *,
    cached_exponent: float,
    minimum: float = 0.0,
    maximum: float = 4.0,
    step: float = _SOURCE_STEP,
) -> NDArray[np.float64]:
    """Evaluate selected VandyCRUISE ``PowerTable.lookup`` semantics.

    The Java table rounds in-range inputs to the nearest precomputed sample
    (ties upward) and evaluates ``Math.pow`` at that sample. NumPy supplies
    the power operation; input quantization preserves the lookup behavior.
    Out-of-range values and exponents different from the cached exponent use
    direct power, as the Java lookup method does.
    """
    array = np.asarray(values, dtype=np.float64)
    q = float(exponent)
    q_cached = float(cached_exponent)
    low = float(minimum)
    high = float(maximum)
    resolution = float(step)
    if not np.isfinite(array).all():
        raise ValueError("PowerTable inputs must be finite")
    if not np.isfinite([q, q_cached, low, high, resolution]).all():
        raise ValueError("PowerTable parameters must be finite")
    if low >= high or resolution <= 0.0:
        raise ValueError("PowerTable range must be ordered and step positive")

    with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
        result = np.power(array, q)
        in_range = (array >= low) & (array <= high)
        if q == q_cached and np.any(in_range):
            indices = np.trunc((array - low) / resolution + 0.5)
            grid = low + indices * resolution
            result = np.where(in_range, np.power(grid, q_cached), result)
    return np.asarray(result, dtype=np.float64)


def toads_optimized_tissue_data_term(
    image_values: ArrayLike,
    centroids: ArrayLike,
    exponent_by_channel_class: ArrayLike,
    modality_weights: ArrayLike,
    class_variance_weights: ArrayLike,
    image_variances: ArrayLike,
    *,
    intensity_field: ArrayLike | None = None,
    psf_bias: ArrayLike | None = None,
) -> NDArray[np.float32]:
    """Compute the selected optimized CSF/GM/WM intensity data terms.

    Inputs have a final modality axis. ``centroids`` and
    ``exponent_by_channel_class`` are ``(channel, class)``; variance weights
    have the same shape, and modality weights/image variances are one value
    per channel. Source float32 residual and weight operations are retained;
    each cached power uses :func:`toads_powertable_lookup`.

    This implements only the source per-class intensity likelihood term. It
    does not add smoothing, atlas priors, lesion/outlier mixtures, membership
    normalization, or topology gates.
    """
    images = np.asarray(image_values, dtype=np.float32)
    means = np.asarray(centroids, dtype=np.float32)
    exponents = np.asarray(exponent_by_channel_class, dtype=np.float32)
    var_weights = np.asarray(class_variance_weights, dtype=np.float32)
    mod_weights = np.asarray(modality_weights, dtype=np.float32)
    var_images = np.asarray(image_variances, dtype=np.float32)
    if images.ndim < 1 or images.shape[-1] < 1 or not np.isfinite(images).all():
        raise ValueError("image_values must be finite with a nonempty modality axis")
    if means.ndim != 2 or means.shape[0] != images.shape[-1] or means.shape[1] < 1:
        raise ValueError("centroids must have shape (modality, class)")
    expected = means.shape
    if exponents.shape != expected or var_weights.shape != expected:
        raise ValueError("exponents and class variance weights must match centroid shape")
    if mod_weights.shape != (images.shape[-1],) or var_images.shape != mod_weights.shape:
        raise ValueError("modality weights and image variances must match modality count")
    if not all(
        np.isfinite(item).all()
        for item in (means, exponents, var_weights, mod_weights, var_images)
    ):
        raise ValueError("TOADS likelihood parameters must be finite")

    if intensity_field is None:
        field_values = np.ones_like(images, dtype=np.float32)
    else:
        field_values = np.asarray(intensity_field, dtype=np.float32)
        if field_values.shape != images.shape or not np.isfinite(field_values).all():
            raise ValueError("intensity_field must be finite and match image_values")
    if psf_bias is None:
        bias_values = np.zeros_like(images, dtype=np.float32)
    else:
        bias_values = np.asarray(psf_bias, dtype=np.float32)
        if bias_values.shape != images.shape or not np.isfinite(bias_values).all():
            raise ValueError("psf_bias must be finite and match image_values")

    corrected = np.asarray(np.asarray(field_values * images, dtype=np.float32) - bias_values, dtype=np.float32)
    energy = np.zeros(images.shape[:-1] + (means.shape[1],), dtype=np.float32)
    for class_index in range(means.shape[1]):
        for channel in range(images.shape[-1]):
            residual = np.asarray(corrected[..., channel] - means[channel, class_index], dtype=np.float32)
            distance = np.asarray(residual * residual, dtype=np.float32)
            powered = toads_powertable_lookup(
                distance,
                float(exponents[channel, class_index]),
                cached_exponent=float(exponents[channel, class_index]),
            )
            weight = np.float32(
                np.float32(mod_weights[channel] * var_weights[channel, class_index])
                * var_images[channel]
            )
            total = energy[..., class_index].astype(np.float64) + powered * np.float64(weight)
            energy[..., class_index] = np.asarray(total, dtype=np.float32)
    return energy
