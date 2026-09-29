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


def toads_ordered_membership_sweep(
    data_energy: ArrayLike,
    initial_memberships: ArrayLike,
    shape_priors: ArrayLike,
    closest_labels: ArrayLike,
    smoothing_priors: ArrayLike,
    smoothing: float,
    prior_factor: float,
    *,
    prior_inverse_energy: ArrayLike | None = None,
    prior_objective_energy: ArrayLike | None = None,
    zero_threshold: float = 1e-10,
    infinite_reciprocal: float = 1e10,
) -> tuple[NDArray[np.float32], np.float32]:
    """Run one source-ordered direct optimized-tissue membership sweep.

    This small-grid stage consumes precomputed per-class data energies. It
    visits interior voxels in x/y/z order and classes in ascending order,
    reads six-face smoothing from the current in-place membership volume,
    adds the half smoothing term to objective energy, applies shape relation
    factors and the source reciprocal/normalization update, and accumulates
    the objective in source order. Optional atlas inverse-prior and objective
    terms are added after smoothing and before reciprocal conversion. It
    omits topology-admissibility checks and does not generate data or PSF terms.
    """
    energies = np.asarray(data_energy, dtype=np.float32)
    initial = np.asarray(initial_memberships, dtype=np.float32)
    smooth_priors = np.asarray(smoothing_priors, dtype=np.float32)
    smooth = np.float32(smoothing)
    threshold = np.float32(zero_threshold)
    reciprocal_cap = np.float32(infinite_reciprocal)
    if energies.ndim != 4 or energies.shape[-1] < 1 or not np.isfinite(energies).all():
        raise ValueError("data_energy must be finite 4-D with a class axis")
    if initial.shape != energies.shape or not np.isfinite(initial).all():
        raise ValueError("initial_memberships must be finite and match data_energy")
    if smooth_priors.shape != (energies.shape[-1],) or not np.isfinite(smooth_priors).all():
        raise ValueError("smoothing_priors must be finite and match the class count")
    if not np.isfinite(smooth) or not np.isfinite(threshold) or not np.isfinite(reciprocal_cap):
        raise ValueError("smoothing and reciprocal constants must be finite")
    if threshold < 0 or reciprocal_cap <= 0:
        raise ValueError("reciprocal threshold/cap must be nonnegative/positive")
    factors, suppressed = toads_shape_relation_factors(shape_priors, closest_labels, prior_factor)
    if factors.shape != energies.shape:
        raise ValueError("shape_priors must match data_energy grid and class count")
    if (prior_inverse_energy is None) != (prior_objective_energy is None):
        raise ValueError("prior inverse and objective energy must be provided together")
    if prior_inverse_energy is None:
        prior_num = None
        prior_objective = None
    else:
        prior_num = np.asarray(prior_inverse_energy, dtype=np.float32)
        prior_objective = np.asarray(prior_objective_energy, dtype=np.float32)
        if prior_num.shape != energies.shape or not np.isfinite(prior_num).all():
            raise ValueError("prior_inverse_energy must be finite and match data_energy")
        if prior_objective.shape != energies.shape or not np.isfinite(prior_objective).all():
            raise ValueError("prior_objective_energy must be finite and match data_energy")

    updated = initial.copy()
    total_objective = np.float32(0.0)
    nx, ny, nz, classes = energies.shape
    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                denominator = np.float32(0.0)
                objective_at_voxel = np.zeros(classes, dtype=np.float32)
                for class_index in range(classes):
                    num = energies[x, y, z, class_index]
                    objective_value = num
                    if smooth > np.float32(0.0):
                        neighbor_sum = np.float32(0.0)
                        for dx, dy, dz in offsets:
                            for other_class in range(classes):
                                if other_class == class_index:
                                    continue
                                neighbor_value = updated[x + dx, y + dy, z + dz, other_class]
                                neighbor_sum = np.float32(
                                    neighbor_sum + np.float32(neighbor_value * neighbor_value)
                                )
                        scaled = np.float32(
                            np.float32(smooth_priors[class_index] * smooth) * neighbor_sum
                        )
                        num = np.float32(num + scaled)
                        energy_increment = np.float32(
                            np.float32(np.float32(np.float32(0.5) * smooth_priors[class_index]) * smooth)
                            * neighbor_sum
                        )
                        objective_value = np.float32(objective_value + energy_increment)

                    if suppressed[x, y, z, class_index]:
                        value = np.float32(0.0)
                    else:
                        if prior_num is not None:
                            num = np.float32(num + prior_num[x, y, z, class_index])
                            objective_value = np.float32(
                                objective_value + prior_objective[x, y, z, class_index]
                            )
                        reciprocal = (
                            np.float32(np.float32(1.0) / num)
                            if num > threshold
                            else reciprocal_cap
                        )
                        value = np.float32(factors[x, y, z, class_index] * reciprocal)
                        objective_at_voxel[class_index] = objective_value
                    updated[x, y, z, class_index] = value
                    denominator = np.float32(denominator + value)

                if denominator > np.float32(0.0):
                    for class_index in range(classes):
                        updated[x, y, z, class_index] = np.float32(
                            updated[x, y, z, class_index] / denominator
                        )
                        contribution = np.float32(
                            np.float32(
                                objective_at_voxel[class_index]
                                * updated[x, y, z, class_index]
                            )
                            * updated[x, y, z, class_index]
                        )
                        total_objective = np.float32(total_objective + contribution)
    return updated, total_objective


def toads_optimized_centroid_update(
    image_values: ArrayLike,
    memberships: ArrayLike,
    centroids: ArrayLike,
    exponents: ArrayLike,
    centroid_smoothness: float,
    *,
    intensity_field: ArrayLike | None = None,
    zero_threshold: float = 1e-10,
) -> NDArray[np.float32]:
    """Update optimized CSF/GM/WM centroids with source damped mode.

    Inputs use ``(x, y, z, modality)`` images/membership spatial grids,
    ``(x, y, z, class)`` memberships, and ``(modality, class)`` centroids and
    exponents. The source accumulates ``factor * membership² * intensity``
    and ``factor * membership²`` in modality/class/x/y/z order, where
    ``factor = PowerTable(distance) / distance`` above ``ZERO`` and one at or
    below it. The result is damped against the prior centroid; a class with
    no positive support keeps its previous centroid. Only interior voxels are
    used, matching the selected source loop.
    """
    images = np.asarray(image_values, dtype=np.float32)
    mems = np.asarray(memberships, dtype=np.float32)
    old_centroids = np.asarray(centroids, dtype=np.float32)
    powers = np.asarray(exponents, dtype=np.float32)
    smooth = np.float32(centroid_smoothness)
    zero = np.float32(zero_threshold)
    if images.ndim != 4 or images.shape[-1] < 1 or not np.isfinite(images).all():
        raise ValueError("image_values must be finite 4-D with a modality axis")
    if mems.ndim != 4 or mems.shape[:3] != images.shape[:3] or mems.shape[-1] < 1:
        raise ValueError("memberships must be 4-D and match image spatial dimensions")
    if not np.isfinite(mems).all():
        raise ValueError("memberships must be finite")
    expected = (images.shape[-1], mems.shape[-1])
    if old_centroids.shape != expected or powers.shape != expected:
        raise ValueError("centroids and exponents must match (modality, class) dimensions")
    if not np.isfinite(old_centroids).all() or not np.isfinite(powers).all():
        raise ValueError("centroids and exponents must be finite")
    if not np.isfinite(smooth) or not np.isfinite(zero) or zero < 0.0:
        raise ValueError("centroid smoothness and ZERO threshold must be finite; ZERO nonnegative")
    if intensity_field is None:
        corrected = images
    else:
        field_values = np.asarray(intensity_field, dtype=np.float32)
        if field_values.shape != images.shape or not np.isfinite(field_values).all():
            raise ValueError("intensity_field must be finite and match image_values")
        corrected = np.asarray(field_values * images, dtype=np.float32)

    updated = old_centroids.copy()
    nx, ny, nz, modalities = images.shape
    classes = mems.shape[-1]
    one_minus_smooth = np.float32(np.float32(1.0) - smooth)
    for channel in range(modalities):
        for class_index in range(classes):
            numerator = np.float32(0.0)
            denominator = np.float32(0.0)
            exponent = float(powers[channel, class_index])
            for x in range(1, nx - 1):
                for y in range(1, ny - 1):
                    for z in range(1, nz - 1):
                        image_value = corrected[x, y, z, channel]
                        residual = np.float32(image_value - old_centroids[channel, class_index])
                        distance = np.float32(residual * residual)
                        factor = np.float32(1.0)
                        if distance > zero:
                            powered = toads_powertable_lookup(
                                distance,
                                exponent,
                                cached_exponent=exponent,
                            )
                            factor = np.float32(np.float64(powered) / np.float64(distance))
                        membership = mems[x, y, z, class_index]
                        weight = np.float32(np.float32(factor * membership) * membership)
                        numerator = np.float32(
                            numerator + np.float32(weight * image_value)
                        )
                        denominator = np.float32(denominator + weight)
            if denominator > np.float32(0.0):
                weighted_mean = np.float32(
                    np.float32(one_minus_smooth * numerator) / denominator
                )
                damped_old = np.float32(smooth * old_centroids[channel, class_index])
                updated[channel, class_index] = np.float32(weighted_mean + damped_old)
    return updated


def toads_optimized_image_variance(
    image_values: ArrayLike,
    memberships: ArrayLike,
    centroids: ArrayLike,
    *,
    intensity_field: ArrayLike | None = None,
    infinite_variance: float = 1e10,
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Compute source image-variance weights for optimized tissue classes.

    Every membership channel is assumed to be an optimized CSF, GM, or WM
    object. The source scans the full grid (including borders), accumulates
    ``(membership**2 + 1e-30) * residual**4`` and
    ``membership**2 + 1e-30`` in class/x/y/z order for each modality, then
    normalizes inverse variances across modalities. The second return value
    is the source's image-level standard deviation, computed from the
    *unnormalized* inverse variance. Optional fields multiply the intensity
    before the residual is calculated.

    Source `INF` behavior is preserved: a channel with no positive numerator
    receives ``infinite_variance`` before cross-modality normalization, so
    all-zero residuals can yield NaN normalized weights as in Java.
    """
    images = np.asarray(image_values, dtype=np.float32)
    mems = np.asarray(memberships, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    if images.ndim != 4 or images.shape[-1] < 1 or not np.isfinite(images).all():
        raise ValueError("image_values must be finite 4-D with a modality axis")
    if mems.ndim != 4 or mems.shape[:3] != images.shape[:3] or mems.shape[-1] < 1:
        raise ValueError("memberships must be finite 4-D and match image spatial dimensions")
    if not np.isfinite(mems).all():
        raise ValueError("memberships must be finite")
    expected = (images.shape[-1], mems.shape[-1])
    if centers.shape != expected or not np.isfinite(centers).all():
        raise ValueError("centroids must be finite and match (modality, class) dimensions")
    infvar = np.float32(infinite_variance)
    if not np.isfinite(infvar) or infvar <= np.float32(0.0):
        raise ValueError("infinite_variance must be finite and positive")
    if intensity_field is None:
        corrected = images
    else:
        field = np.asarray(intensity_field, dtype=np.float32)
        if field.shape != images.shape or not np.isfinite(field).all():
            raise ValueError("intensity_field must be finite and match image_values")
        corrected = np.asarray(field * images, dtype=np.float32)

    modalities = images.shape[-1]
    classes = mems.shape[-1]
    variances = np.empty(modalities, dtype=np.float32)
    for channel in range(modalities):
        numerator = np.float32(0.0)
        denominator = np.float32(0.0)
        for class_index in range(classes):
            for x in range(images.shape[0]):
                for y in range(images.shape[1]):
                    for z in range(images.shape[2]):
                        membership = mems[x, y, z, class_index]
                        weight = np.float32(
                            np.float32(membership * membership) + np.float32(1e-30)
                        )
                        residual = np.float32(
                            corrected[x, y, z, channel] - centers[channel, class_index]
                        )
                        difference = np.float32(residual * residual)
                        weighted_difference = np.float32(weight * difference)
                        contribution = np.float32(weighted_difference * difference)
                        numerator = np.float32(
                            numerator + contribution
                        )
                        denominator = np.float32(denominator + weight)
        if numerator > np.float32(0.0):
            variances[channel] = np.float32(denominator / numerator)
        else:
            variances[channel] = infvar

    total = np.float32(0.0)
    for value in variances:
        total = np.float32(total + value)
    image_stddev = np.asarray(np.sqrt(variances), dtype=np.float32)
    with np.errstate(divide="ignore", invalid="ignore"):
        normalized = np.asarray(variances / total, dtype=np.float32)
    return normalized, image_stddev


def toads_optimized_class_standard_deviation(
    image_values: ArrayLike,
    memberships: ArrayLike,
    centroids: ArrayLike,
    segmentation: ArrayLike,
    *,
    intensity_field: ArrayLike | None = None,
) -> NDArray[np.float32]:
    """Compute assigned-class standard deviations for optimized TOADS types.

    This transcribes the optimized tissue branch of source
    `computeExactVariances`. Each modality/class uses only voxels whose
    categorical segmentation equals that class. The source includes the full
    image grid (borders included), uses ``membership**2 + 1e-30`` as its
    weight, and returns zero when no voxel is assigned to a class. The result
    has shape ``(modality, class)``. All supplied classes are assumed to be
    optimized CSF, GM, or WM types.
    """
    images = np.asarray(image_values, dtype=np.float32)
    mems = np.asarray(memberships, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    labels = np.asarray(segmentation)
    if images.ndim != 4 or images.shape[-1] < 1 or not np.isfinite(images).all():
        raise ValueError("image_values must be finite 4-D with a modality axis")
    if mems.ndim != 4 or mems.shape[:3] != images.shape[:3] or mems.shape[-1] < 1:
        raise ValueError("memberships must be finite 4-D and match image spatial dimensions")
    if not np.isfinite(mems).all():
        raise ValueError("memberships must be finite")
    if labels.shape != images.shape[:3] or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer array matching the image grid")
    expected = (images.shape[-1], mems.shape[-1])
    if centers.shape != expected or not np.isfinite(centers).all():
        raise ValueError("centroids must be finite and match (modality, class) dimensions")
    if intensity_field is None:
        corrected = images
    else:
        field = np.asarray(intensity_field, dtype=np.float32)
        if field.shape != images.shape or not np.isfinite(field).all():
            raise ValueError("intensity_field must be finite and match image_values")
        corrected = np.asarray(field * images, dtype=np.float32)

    deviations = np.zeros(expected, dtype=np.float32)
    for channel in range(images.shape[-1]):
        for class_index in range(mems.shape[-1]):
            numerator = np.float32(0.0)
            denominator = np.float32(0.0)
            for x in range(images.shape[0]):
                for y in range(images.shape[1]):
                    for z in range(images.shape[2]):
                        if labels[x, y, z] != class_index:
                            continue
                        membership = mems[x, y, z, class_index]
                        membership_squared = np.float32(membership * membership)
                        weight = np.float32(membership_squared + np.float32(1e-30))
                        residual = np.float32(
                            corrected[x, y, z, channel] - centers[channel, class_index]
                        )
                        weighted_residual = np.float32(weight * residual)
                        contribution = np.float32(weighted_residual * residual)
                        numerator = np.float32(
                            numerator + contribution
                        )
                        denominator = np.float32(denominator + weight)
            if denominator > np.float32(0.0):
                deviations[channel, class_index] = np.sqrt(
                    np.float32(numerator / denominator)
                )
    return deviations


def toads_optimized_centroid_similarity(
    centroids: ArrayLike,
    class_standard_deviations: ArrayLike,
    image_variance_weights: ArrayLike,
    prior_scale: float = 0.05,
) -> NDArray[np.float32]:
    """Compute source centroid-similarity coefficients for optimized classes.

    The source accumulates a modality-weighted squared centroid distance,
    divided twice by the sum of the two class standard deviations, then maps
    it through ``1 / (1 + (1-priorScale)/priorScale * distance)``. Diagonal
    entries are explicitly zero. Float32 operation and modality order match
    the selected `computeExactMemberships` formula. IEEE nonfinite results
    from zero deviation sums or singular ``prior_scale`` values are preserved.
    This is the pairwise coefficient only; it does not apply voxelwise shape
    priors.
    """
    centers = np.asarray(centroids, dtype=np.float32)
    deviations = np.asarray(class_standard_deviations, dtype=np.float32)
    weights = np.asarray(image_variance_weights, dtype=np.float32)
    scale = np.float32(prior_scale)
    if centers.ndim != 2 or centers.shape[0] < 1 or centers.shape[1] < 1:
        raise ValueError("centroids must be a nonempty (modality, class) array")
    if deviations.shape != centers.shape or not np.isfinite(deviations).all():
        raise ValueError("class_standard_deviations must be finite and match centroids")
    if weights.shape != (centers.shape[0],) or not np.isfinite(weights).all():
        raise ValueError("image_variance_weights must be finite and match modality count")
    if not np.isfinite(centers).all() or not np.isfinite(scale):
        raise ValueError("centroids and prior_scale must be finite")

    classes = centers.shape[1]
    similarities = np.zeros((classes, classes), dtype=np.float32)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        scale_factor = np.float32(np.float32(np.float32(1.0) - scale) / scale)
        for k in range(classes):
            for m in range(classes):
                if m == k:
                    continue
                distance = np.float32(0.0)
                for channel in range(centers.shape[0]):
                    centroid_delta = np.float32(centers[channel, k] - centers[channel, m])
                    deviation_sum = np.float32(
                        deviations[channel, k] + deviations[channel, m]
                    )
                    term = np.float32(weights[channel] * centroid_delta)
                    term = np.float32(term * centroid_delta)
                    term = np.float32(term / deviation_sum)
                    term = np.float32(term / deviation_sum)
                    distance = np.float32(distance + term)
                scaled_distance = np.float32(scale_factor * distance)
                denominator = np.float32(np.float32(1.0) + scaled_distance)
                similarities[k, m] = np.float32(np.float32(1.0) / denominator)
    return similarities


def toads_inverse_shape_prior_energy(
    shape_values: ArrayLike,
    shape_priors: ArrayLike,
    centroid_similarity: ArrayLike,
    prior_coefficient: float,
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Compute TOADS' atlas inverse-prior and objective-energy terms.

    For each interior voxel and candidate class ``k``, the source sums over
    all ``m != k`` the product ``shape_prior[m] * similarity[k,m] *
    shape[m] * shape[m]``. It adds ``prior_coefficient * invprior`` to the
    class inverse-energy and half that amount to objective energy. The source
    only enters this calculation when the coefficient is positive; borders
    therefore receive zero terms here. Class and spatial traversal use the
    selected source order and float32 arithmetic.
    """
    shapes = np.asarray(shape_values, dtype=np.float32)
    priors = np.asarray(shape_priors, dtype=np.float32)
    similarity = np.asarray(centroid_similarity, dtype=np.float32)
    coefficient = np.float32(prior_coefficient)
    if shapes.ndim != 4 or shapes.shape[-1] < 1 or not np.isfinite(shapes).all():
        raise ValueError("shape_values must be finite 4-D with a class axis")
    classes = shapes.shape[-1]
    if priors.shape != (classes,) or not np.isfinite(priors).all():
        raise ValueError("shape_priors must be finite and match the class count")
    if similarity.shape != (classes, classes) or not np.isfinite(similarity).all():
        raise ValueError("centroid_similarity must be finite and square by class")
    if not np.isfinite(coefficient):
        raise ValueError("prior_coefficient must be finite")
    inverse_energy = np.zeros(shapes.shape, dtype=np.float32)
    objective_energy = np.zeros(shapes.shape, dtype=np.float32)
    if coefficient <= np.float32(0.0):
        return inverse_energy, objective_energy

    nx, ny, nz, _ = shapes.shape
    half_coefficient = np.float32(np.float32(0.5) * coefficient)
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                for candidate in range(classes):
                    inverse_prior = np.float32(0.0)
                    for other in range(classes):
                        if other == candidate:
                            continue
                        term = np.float32(priors[other] * similarity[candidate, other])
                        term = np.float32(term * shapes[x, y, z, other])
                        term = np.float32(term * shapes[x, y, z, other])
                        inverse_prior = np.float32(inverse_prior + term)
                    inverse_energy[x, y, z, candidate] = np.float32(
                        coefficient * inverse_prior
                    )
                    objective_energy[x, y, z, candidate] = np.float32(
                        half_coefficient * inverse_prior
                    )
    return inverse_energy, objective_energy
