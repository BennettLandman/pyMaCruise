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


def toads_psf_ridge_bias(
    image_values: ArrayLike,
    segmentation: ArrayLike,
    psf_ridge: float,
    *,
    intensity_field: ArrayLike | None = None,
) -> NDArray[np.float32]:
    """Compute source `computeExactMemberships` PSF bias on 3-D boundaries.

    The source only evaluates interior voxels. A bias is produced where any
    segmentation value in the 3×3×3 neighborhood differs from the center;
    the center sample is weighted by ``-27/26`` and all 27 samples (including
    the center) are added at ``1/26`` before multiplication by ``psf_ridge``.
    Source float32 operation and accumulation order are retained. The optional
    field multiplies each image sample before the stencil.
    """
    images = np.asarray(image_values, dtype=np.float32)
    labels = np.asarray(segmentation)
    if images.ndim != 4 or images.shape[-1] < 1 or not np.isfinite(images).all():
        raise ValueError("image_values must be finite 4-D with a modality axis")
    if labels.shape != images.shape[:3] or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer array matching the image grid")
    ridge = np.float32(psf_ridge)
    if not np.isfinite(ridge):
        raise ValueError("psf_ridge must be finite")
    if intensity_field is None:
        corrected = images
    else:
        field_values = np.asarray(intensity_field, dtype=np.float32)
        if field_values.shape != images.shape or not np.isfinite(field_values).all():
            raise ValueError("intensity_field must be finite and match image_values")
        corrected = np.asarray(field_values * images, dtype=np.float32)

    bias = np.zeros(images.shape, dtype=np.float32)
    if ridge <= np.float32(0.0):
        return bias
    nx, ny, nz, channels = images.shape
    center_weight = np.float32(np.float32(-27.0) / np.float32(26.0))
    neighbor_weight = np.float32(np.float32(1.0) / np.float32(26.0))
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                center_label = labels[x, y, z]
                patch = labels[x - 1 : x + 2, y - 1 : y + 2, z - 1 : z + 2]
                if not np.any(patch != center_label):
                    continue
                for channel in range(channels):
                    value = np.float32(center_weight * corrected[x, y, z, channel])
                    for i in range(-1, 2):
                        for j in range(-1, 2):
                            for k in range(-1, 2):
                                contribution = np.float32(
                                    neighbor_weight * corrected[x + i, y + j, z + k, channel]
                                )
                                value = np.float32(value + contribution)
                    bias[x, y, z, channel] = np.float32(value * ridge)
    return bias


def toads_six_neighbor_smoothing_terms(
    memberships: ArrayLike,
    smoothing: float,
    smoothing_priors: ArrayLike,
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Compute the selected TOADS six-face membership-smoothing terms.

    The input has a final class axis. For each interior voxel and candidate
    class, the source visits X+, X-, Y+, Y-, Z+, Z- neighbors in that order,
    summing squared memberships for every class except the candidate. The
    returned arrays are the increment to reciprocal-energy ``num`` and the
    half-scaled increment to the objective energy. Border voxels are not
    visited by the source update and therefore receive zero terms. The
    source's ``smoothing > 0`` guard also makes nonpositive values return zero.
    """
    mems = np.asarray(memberships, dtype=np.float32)
    priors = np.asarray(smoothing_priors, dtype=np.float32)
    smooth = np.float32(smoothing)
    if mems.ndim != 4 or mems.shape[-1] < 1 or not np.isfinite(mems).all():
        raise ValueError("memberships must be finite 4-D with a class axis")
    if priors.shape != (mems.shape[-1],) or not np.isfinite(priors).all():
        raise ValueError("smoothing_priors must be finite and match the class count")
    if not np.isfinite(smooth):
        raise ValueError("smoothing must be finite")

    num_increment = np.zeros(mems.shape, dtype=np.float32)
    energy_increment = np.zeros(mems.shape, dtype=np.float32)
    if smooth <= np.float32(0.0):
        return num_increment, energy_increment
    nx, ny, nz, classes = mems.shape
    neighbor_offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                for candidate in range(classes):
                    neighbor_sum = np.float32(0.0)
                    for dx, dy, dz in neighbor_offsets:
                        for other_class in range(classes):
                            if other_class == candidate:
                                continue
                            value = mems[x + dx, y + dy, z + dz, other_class]
                            neighbor_sum = np.float32(
                                neighbor_sum + np.float32(value * value)
                            )
                    scaled = np.float32(np.float32(priors[candidate] * smooth) * neighbor_sum)
                    num_increment[x, y, z, candidate] = scaled
                    energy_increment[x, y, z, candidate] = np.float32(
                        np.float32(np.float32(np.float32(0.5) * priors[candidate]) * smooth)
                        * neighbor_sum
                    )
    return num_increment, energy_increment


def toads_direct_membership_update(
    data_energy: ArrayLike,
    initial_memberships: ArrayLike,
    *,
    objective_energy: ArrayLike | None = None,
    zero_threshold: float = 1e-10,
    infinite_reciprocal: float = 1e10,
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Normalize reciprocal energies for direct optimized-tissue updates.

    This captures the selected source's direct relation update after all
    non-data terms have already been added to ``data_energy``. Interior
    values above ``zero_threshold`` are inverted; values at or below it use
    ``infinite_reciprocal``. Class reciprocals are summed and normalized in
    source order. Border memberships are copied unchanged. The second return
    value contains each class's contribution to the source objective sum
    ``objective_energy * membership * membership`` on updated interior
    voxels; border contributions remain zero.

    This is not the asynchronous spatial-smoothing/update sweep or its
    topology/prior relations; those require separate source-stage composition.
    """
    energies = np.asarray(data_energy, dtype=np.float32)
    initial = np.asarray(initial_memberships, dtype=np.float32)
    objective = energies if objective_energy is None else np.asarray(objective_energy, dtype=np.float32)
    threshold = np.float32(zero_threshold)
    reciprocal_cap = np.float32(infinite_reciprocal)
    if energies.ndim != 4 or energies.shape[-1] < 1 or not np.isfinite(energies).all():
        raise ValueError("data_energy must be finite 4-D with a class axis")
    if initial.shape != energies.shape or not np.isfinite(initial).all():
        raise ValueError("initial_memberships must be finite and match data_energy")
    if objective.shape != energies.shape or not np.isfinite(objective).all():
        raise ValueError("objective_energy must be finite and match data_energy")
    if not np.isfinite(threshold) or not np.isfinite(reciprocal_cap) or threshold < 0 or reciprocal_cap <= 0:
        raise ValueError("reciprocal thresholds must be finite and nonnegative/positive")

    updated = initial.copy()
    contributions = np.zeros(energies.shape, dtype=np.float32)
    nx, ny, nz, classes = energies.shape
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                denominator = np.float32(0.0)
                for class_index in range(classes):
                    value = energies[x, y, z, class_index]
                    reciprocal = (
                        np.float32(np.float32(1.0) / value)
                        if value > threshold
                        else reciprocal_cap
                    )
                    updated[x, y, z, class_index] = reciprocal
                    denominator = np.float32(denominator + reciprocal)
                if denominator > np.float32(0.0):
                    for class_index in range(classes):
                        updated[x, y, z, class_index] = np.float32(
                            updated[x, y, z, class_index] / denominator
                        )
                        contribution = np.float32(
                            np.float32(
                                objective[x, y, z, class_index]
                                * updated[x, y, z, class_index]
                            )
                            * updated[x, y, z, class_index]
                        )
                        contributions[x, y, z, class_index] = contribution
    return updated, contributions


def toads_shape_relation_factors(
    shape_priors: ArrayLike,
    closest_labels: ArrayLike,
    prior_factor: float,
) -> tuple[NDArray[np.float32], NDArray[np.bool_]]:
    """Return source TOADS direct/neighbor/prior membership factors.

    ``closest_labels[..., :2]`` are the source direct slots and
    ``closest_labels[..., 2:]`` are its neighbor slots. Direct relation has
    precedence. Neighbor factors use the source clamp
    ``max(0.25, min(2*prior_factor*shape, 1))``; prior factors use
    ``min(1, 2*prior_factor*shape)``. A prior-only class with shape exactly
    zero is separately marked as suppressed, matching the update's early
    support check.
    """
    shapes = np.asarray(shape_priors, dtype=np.float32)
    labels = np.asarray(closest_labels)
    factor = np.float32(prior_factor)
    if shapes.ndim != 4 or shapes.shape[-1] < 1 or not np.isfinite(shapes).all():
        raise ValueError("shape_priors must be finite 4-D with a class axis")
    if labels.shape != shapes.shape[:3] + (4,) or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("closest_labels must be an integer four-slot array matching the shape grid")
    if not np.isfinite(factor):
        raise ValueError("prior_factor must be finite")
    classes = shapes.shape[-1]
    if np.any((labels < -1) | (labels >= classes)):
        raise ValueError("closest_labels values must be -1 or valid class indices")

    factors = np.empty(shapes.shape, dtype=np.float32)
    suppressed = np.zeros(shapes.shape, dtype=np.bool_)
    for x in range(shapes.shape[0]):
        for y in range(shapes.shape[1]):
            for z in range(shapes.shape[2]):
                for class_index in range(classes):
                    shape = shapes[x, y, z, class_index]
                    direct = bool(np.any(labels[x, y, z, :2] == class_index))
                    neighbor = bool(np.any(labels[x, y, z, 2:] == class_index))
                    if direct:
                        factors[x, y, z, class_index] = np.float32(1.0)
                    elif neighbor:
                        scaled = np.float32(np.float32(2.0) * factor)
                        scaled = np.float32(scaled * shape)
                        factors[x, y, z, class_index] = np.float32(
                            max(np.float32(0.25), min(scaled, np.float32(1.0)))
                        )
                    elif shape == np.float32(0.0):
                        factors[x, y, z, class_index] = np.float32(0.0)
                        suppressed[x, y, z, class_index] = True
                    else:
                        scaled = np.float32(np.float32(2.0) * factor)
                        scaled = np.float32(scaled * shape)
                        factors[x, y, z, class_index] = np.float32(
                            min(np.float32(1.0), scaled)
                        )
    return factors, suppressed


def toads_shape_weighted_membership_update(
    data_energy: ArrayLike,
    initial_memberships: ArrayLike,
    shape_priors: ArrayLike,
    closest_labels: ArrayLike,
    prior_factor: float,
    *,
    zero_threshold: float = 1e-10,
    infinite_reciprocal: float = 1e10,
) -> NDArray[np.float32]:
    """Apply source relation factors to reciprocals, then normalize classes.

    This composes the direct/neighbor/prior relation gate with the selected
    optimized-tissue reciprocal-energy update at each interior voxel. A
    prior-only class with zero shape is set to zero before the denominator is
    summed. Smoothing and other terms are expected to have been composed into
    ``data_energy`` by a caller; this helper does not implement the source's
    asynchronous volume sweep.
    """
    energies = np.asarray(data_energy, dtype=np.float32)
    initial = np.asarray(initial_memberships, dtype=np.float32)
    if energies.ndim != 4 or energies.shape[-1] < 1 or not np.isfinite(energies).all():
        raise ValueError("data_energy must be finite 4-D with a class axis")
    if initial.shape != energies.shape or not np.isfinite(initial).all():
        raise ValueError("initial_memberships must be finite and match data_energy")
    factors, suppressed = toads_shape_relation_factors(shape_priors, closest_labels, prior_factor)
    if factors.shape != energies.shape:
        raise ValueError("shape_priors must match data_energy grid and class count")
    threshold = np.float32(zero_threshold)
    reciprocal_cap = np.float32(infinite_reciprocal)
    if not np.isfinite(threshold) or not np.isfinite(reciprocal_cap) or threshold < 0 or reciprocal_cap <= 0:
        raise ValueError("reciprocal thresholds must be finite and nonnegative/positive")

    updated = initial.copy()
    nx, ny, nz, classes = energies.shape
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                denominator = np.float32(0.0)
                for class_index in range(classes):
                    if suppressed[x, y, z, class_index]:
                        value = np.float32(0.0)
                    else:
                        energy = energies[x, y, z, class_index]
                        reciprocal = (
                            np.float32(np.float32(1.0) / energy)
                            if energy > threshold
                            else reciprocal_cap
                        )
                        value = np.float32(factors[x, y, z, class_index] * reciprocal)
                    updated[x, y, z, class_index] = value
                    denominator = np.float32(denominator + value)
                if denominator > np.float32(0.0):
                    for class_index in range(classes):
                        updated[x, y, z, class_index] = np.float32(
                            updated[x, y, z, class_index] / denominator
                        )
    return updated
