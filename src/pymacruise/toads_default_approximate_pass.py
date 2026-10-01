"""Alternate-constructor approx=4 update for the default TOADS atlas.

This follows VandyCRUISE ``MultipleToadSegmentation.computeApproxMemberships``
for the default 11-class atlas when called through the constructor that sets
``approx=4``. The shorter constructor used by ``MedicAlgorithmToads08`` forces
``useApprox=false``; this API is not evidence that that wrapper runs this mode.
The source's top-k sparse state and in-place x/y/z sweep are retained; this is
not a full TOADS loop.

Paper context: Bazin and Pham, “Homeomorphic brain image segmentation with
topological and statistical atlases” (2008, DOI 10.1016/j.media.2008.06.008),
§2.1, Eqs. 1–4 describe fuzzy tissue memberships with gain-field-corrected
intensity terms, neighborhood smoothing, atlas/relation contributions, and
class intensity deviations used to shape atlas influence. This pass has
conceptual correspondence through its approximate image-variance and
centroid-similarity preludes and its data/smoothing/prior score terms. Its
four-slot approximation, selected score equations and dispatch,
candidate normalization, constructor profile, and objective details follow
VandyCRUISE source; they are not claims about the published equations.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES
from .toads_math import toads_powertable_lookup

_ZERO = np.float32(1.0e-10)
_INF = np.float32(1.0e10)
_EPS = np.float32(1.0e-30)


@dataclass(frozen=True)
class ToadsDefaultApproximatePass:
    """Updated sparse state and method-level diagnostics."""

    memberships: NDArray[np.float32]
    best_class_indices: NDArray[np.uint8]
    denominator: NDArray[np.float32]
    image_variances: NDArray[np.float32]
    image_variance_numerator: NDArray[np.float32]
    image_variance_denominator: NDArray[np.float32]
    unnormalized_image_variance: NDArray[np.float32]
    image_standard_deviation: NDArray[np.float32]
    centroid_similarity: NDArray[np.float32]
    objective: np.float32


def _f(value: float | np.floating) -> np.float32:
    return np.float32(value)


def _similarity(
    centers: NDArray[np.float32],
    deviations: NDArray[np.float32],
    image_variance: NDArray[np.float32],
    prior_scale: np.float32,
) -> NDArray[np.float32]:
    channels, classes = centers.shape
    result = np.zeros((classes, classes), dtype=np.float32)
    for k in range(classes):
        for m in range(classes):
            if m == k:
                continue
            distance = _f(0.0)
            for c in range(channels):
                delta = _f(centers[c, k] - centers[c, m])
                term = _f(image_variance[c] * delta)
                term = _f(term * delta)
                denominator = _f(deviations[c, k] + deviations[c, m])
                term = _f(term / denominator)
                term = _f(term / denominator)
                distance = _f(distance + term)
            scale = _f(_f(_f(1.0) - prior_scale) / prior_scale)
            result[k, m] = _f(_f(1.0) / _f(_f(1.0) + _f(scale * distance)))
    return result


def _approximate_image_variance(
    images: NDArray[np.float32],
    memberships: NDArray[np.float32],
    best: NDArray[np.uint8],
    centers: NDArray[np.float32],
    field: NDArray[np.float32] | None,
) -> tuple[NDArray[np.float32], NDArray[np.float32], NDArray[np.float32], NDArray[np.float32], NDArray[np.float32]]:
    nx, ny, nz, channels = images.shape
    classes = centers.shape[1]
    num = np.zeros(channels, dtype=np.float32)
    den = np.zeros(channels, dtype=np.float32)
    for c in range(channels):
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    for a in range(4):
                        k = int(best[x, y, z, a])
                        # The selected default atlas has MASK at index 0;
                        # only optimized image classes contribute here.
                        if k == 0 or not 0 <= k < classes:
                            continue
                        value = images[x, y, z, c]
                        if field is not None:
                            value = _f(field[x, y, z, c] * value)
                        residual = _f(value - centers[c, k])
                        squared = _f(residual * residual)
                        weight = _f(_f(memberships[x, y, z, a] * memberships[x, y, z, a]) + _EPS)
                        contribution = _f(weight * squared)
                        contribution = _f(contribution * squared)
                        num[c] = _f(num[c] + contribution)
                        den[c] = _f(den[c] + weight)
    raw = np.empty(channels, dtype=np.float32)
    total = _f(0.0)
    for c in range(channels):
        raw[c] = _f(den[c] / num[c]) if num[c] > _f(0.0) else _f(np.inf)
        total = _f(total + raw[c])
    output = np.empty(channels, dtype=np.float32)
    deviation = np.empty(channels, dtype=np.float32)
    for c in range(channels):
        output[c] = _f(raw[c] / total)
        deviation[c] = _f(np.sqrt(np.float64(raw[c])))
    return output, num, den, raw, deviation


def toads_default_approximate_membership_pass(
    image_values: ArrayLike,
    initial_memberships: ArrayLike,
    best_class_indices: ArrayLike,
    centroids: ArrayLike,
    standard_deviations: ArrayLike,
    class_variance_weights: ArrayLike,
    modality_weights: ArrayLike,
    exponent_by_class: ArrayLike,
    transformed_shapes: ArrayLike,
    direct_relation: ArrayLike,
    smoothing_priors: ArrayLike,
    atlas_shape_priors: ArrayLike,
    *,
    prior_scale: float,
    prior_factor: float,
    prior_coefficient: float,
    smoothing: float,
    intensity_field: ArrayLike | None = None,
) -> ToadsDefaultApproximatePass:
    """Run the alternate four-slot approximate membership update.

    This alternate constructor path uses four candidates with the default
    11-class atlas. Fixed upstream inputs are in true-RAS voxel order. ``direct_relation``
    represents ``relationList[relationMap[x][y][z]][k]``; transformed atlas
    shapes and statistics are supplied by prior stages. The method recomputes
    the selected approximate image-variance and centroid-similarity preludes.
    """
    images = np.asarray(image_values, dtype=np.float32)
    memberships = np.asarray(initial_memberships, dtype=np.float32).copy()
    best_raw = np.asarray(best_class_indices)
    centers = np.asarray(centroids, dtype=np.float32)
    deviations = np.asarray(standard_deviations, dtype=np.float32)
    variance_weights = np.asarray(class_variance_weights, dtype=np.float32)
    modality = np.asarray(modality_weights, dtype=np.float32)
    exponents = np.asarray(exponent_by_class, dtype=np.float32)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    direct = np.asarray(direct_relation, dtype=bool)
    smooth_prior = np.asarray(smoothing_priors, dtype=np.float32)
    shape_prior = np.asarray(atlas_shape_priors, dtype=np.float32)
    types = DEFAULT_ATLAS_CLASS_TYPES
    nx, ny, nz, channels = images.shape if images.ndim == 4 else (0, 0, 0, 0)
    spatial, classes = (nx, ny, nz), len(types)
    if channels != 1 or min(spatial) < 3 or memberships.shape != spatial + (4,):
        raise ValueError("selected default pass requires (x>=3,y>=3,z>=3,1) and four slots")
    if best_raw.shape != memberships.shape or not np.issubdtype(best_raw.dtype, np.integer):
        raise ValueError("best_class_indices must be an integer four-slot grid")
    if centers.shape != (channels, classes) or deviations.shape != centers.shape:
        raise ValueError("centroids and standard_deviations must have shape (1,11)")
    if variance_weights.shape != centers.shape or exponents.shape != (classes,):
        raise ValueError("class weights must be (1,11) and exponents must have 11 values")
    if modality.shape != (channels,) or shapes.shape != spatial + (classes,):
        raise ValueError("modality weights and transformed shapes have incompatible dimensions")
    if direct.shape != shapes.shape or smooth_prior.shape != (classes,) or shape_prior.shape != (classes,):
        raise ValueError("relation and atlas-prior arrays must match the 11-class grid")
    arrays = (images, memberships, centers, deviations, variance_weights, modality,
              exponents, shapes, smooth_prior, shape_prior)
    if not all(np.isfinite(value).all() for value in arrays):
        raise ValueError("fixed TOADS inputs must be finite")
    if np.any((best_raw < 0) | (best_raw >= classes)):
        raise ValueError("candidate class indices must be in [0,10]")
    best = best_raw.astype(np.uint8, copy=True)
    field = None if intensity_field is None else np.asarray(intensity_field, dtype=np.float32)
    if field is not None and (field.shape != images.shape or not np.isfinite(field).all()):
        raise ValueError("intensity_field must be finite and match the one-channel image")
    prior_scale_f = _f(prior_scale)
    prior_factor_f = _f(prior_factor)
    prior_coefficient_f = _f(prior_coefficient)
    smoothing_f = _f(smoothing)
    if not np.isfinite([prior_scale_f, prior_factor_f, prior_coefficient_f, smoothing_f]).all():
        raise ValueError("TOADS scalar settings must be finite")
    if prior_scale_f == 0 or smoothing_f < 0:
        raise ValueError("prior_scale must be nonzero and smoothing nonnegative")

    image_variance, image_variance_num, image_variance_den, raw_image_variance, image_std = _approximate_image_variance(
        images, memberships, best, centers, field
    )
    similarity = _similarity(centers, deviations, image_variance, prior_scale_f)
    denominator_map = np.zeros(spatial, dtype=np.float32)
    objective = _f(0.0)
    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                voxel_shape = shapes[x, y, z] if prior_coefficient_f > _f(0.0) else np.zeros(classes, np.float32)
                energy = np.zeros(classes, dtype=np.float32)
                scores = np.zeros(classes, dtype=np.float32)
                denominator = _f(0.0)
                for k, kind in enumerate(types):
                    if not direct[x, y, z, k] and voxel_shape[k] == _f(0.0):
                        continue
                    image = images[x, y, z, 0]
                    if kind == "MASK":
                        residual_cost = _f(_f(modality[0] * variance_weights[0, k]) * image_variance[0])
                        if image != _f(0.0):
                            residual_cost = _f(residual_cost * _f(1.0))
                        else:
                            residual_cost = _f(0.0)
                    else:
                        corrected = _f(image if field is None else _f(field[x, y, z, 0] * image))
                        residual = _f(corrected - centers[0, k])
                        distance = _f(residual * residual)
                        powered = toads_powertable_lookup(
                            np.asarray(distance, dtype=np.float64),
                            float(exponents[k]), cached_exponent=float(exponents[k]),
                        ).item()
                        weight = _f(_f(modality[0] * variance_weights[0, k]) * image_variance[0])
                        residual_cost = _f(float(weight) * float(powered))
                    energy[k] = residual_cost
                    num = residual_cost
                    if smoothing_f > _f(0.0):
                        neighbor_sum = _f(0.0)
                        for dx, dy, dz in offsets:
                            qx, qy, qz = x + dx, y + dy, z + dz
                            for a in range(4):
                                if int(best[qx, qy, qz, a]) != k:
                                    value = memberships[qx, qy, qz, a]
                                    neighbor_sum = _f(neighbor_sum + _f(value * value))
                        increment = _f(_f(_f(smooth_prior[k] * smoothing_f) * neighbor_sum))
                        num = _f(num + increment)
                        energy[k] = _f(energy[k] + _f(_f(_f(_f(0.5) * smooth_prior[k]) * smoothing_f) * neighbor_sum))
                    if prior_coefficient_f > _f(0.0):
                        inverse_prior = _f(0.0)
                        for m in range(classes):
                            if m == k:
                                continue
                            term = _f(shape_prior[m] * similarity[k, m])
                            term = _f(term * voxel_shape[m])
                            term = _f(term * voxel_shape[m])
                            inverse_prior = _f(inverse_prior + term)
                        num = _f(num + _f(prior_coefficient_f * inverse_prior))
                        energy[k] = _f(energy[k] + _f(_f(_f(0.5) * prior_coefficient_f) * inverse_prior))
                    score = _f(_f(1.0) / num) if num > _ZERO else _INF
                    if not direct[x, y, z, k]:
                        shape_factor = min(_f(1.0), _f(_f(_f(2.0) * prior_factor_f) * voxel_shape[k]))
                        score = _f(shape_factor * score)
                    scores[k] = score
                    denominator = _f(denominator + score)
                denominator_map[x, y, z] = denominator
                scratch = scores.copy()
                for a in range(4):
                    winner = 0
                    for candidate in range(1, classes):
                        if scratch[candidate] > scratch[winner]:
                            winner = candidate
                    best[x, y, z, a] = winner
                    memberships[x, y, z, a] = scratch[winner]
                    scratch[winner] = -_INF
                if denominator > _f(0.0):
                    for a in range(4):
                        k = int(best[x, y, z, a])
                        memberships[x, y, z, a] = _f(memberships[x, y, z, a] / denominator)
                        objective = _f(objective + _f(_f(energy[k] * memberships[x, y, z, a]) * memberships[x, y, z, a]))
    return ToadsDefaultApproximatePass(
        memberships, best, denominator_map, image_variance,
        image_variance_num, image_variance_den, raw_image_variance, image_std,
        similarity, objective,
    )
