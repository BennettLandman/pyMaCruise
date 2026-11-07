"""Selected fixed-input AtlasEM statistics and membership updates.

This module reproduces the active ``AtlasEMSegmentation.computePriors``,
``computeCentriods`` (source spelling), and ``computeCovariance`` loops. It is
not the AtlasEM pipeline: full preprocessing, atlas registration, complete EM
iteration, and the missing selected atlas files are tracked under P44b.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from numba import njit


@njit(cache=True, fastmath=False)
def _centroids_kernel(values, mem, active, field, use_field, means):
    for c in range(values.shape[3]):
        for k in range(mem.shape[3]):
            numerator = np.float32(0.0)
            denominator = np.float32(0.0)
            for x in range(values.shape[0]):
                for y in range(values.shape[1]):
                    for z in range(values.shape[2]):
                        if active[x, y, z]:
                            membership = mem[x, y, z, k]
                            sample = values[x, y, z, c]
                            if use_field:
                                weighted = np.float32(membership * field[x, y, z, c])
                                weighted = np.float32(weighted * sample)
                            else:
                                weighted = np.float32(membership * sample)
                            numerator = np.float32(numerator + weighted)
                            denominator = np.float32(denominator + membership)
            if denominator != np.float32(0.0):
                means[c, k] = np.float32(numerator / denominator)
    return means


@njit(cache=True, fastmath=False)
def _covariances_kernel(values, mem, active, means, field, use_field):
    channels, classes = values.shape[3], mem.shape[3]
    covariance = np.zeros((classes, channels, channels), dtype=np.float32)
    for k in range(classes):
        denominator = np.float32(0.0)
        for x in range(values.shape[0]):
            for y in range(values.shape[1]):
                for z in range(values.shape[2]):
                    if active[x, y, z]:
                        membership = mem[x, y, z, k]
                        for c1 in range(channels):
                            sample1 = values[x, y, z, c1]
                            if use_field:
                                sample1 = np.float32(field[x, y, z, c1] * sample1)
                            d1 = np.float32(sample1 - means[c1, k])
                            for c2 in range(channels):
                                sample2 = values[x, y, z, c2]
                                if use_field:
                                    sample2 = np.float32(field[x, y, z, c2] * sample2)
                                d2 = np.float32(sample2 - means[c2, k])
                                product = np.float32(np.float32(d1 * d2) * membership)
                                covariance[k, c1, c2] = np.float32(covariance[k, c1, c2] + product)
                                denominator = np.float32(denominator + membership)
        if denominator != np.float32(0.0):
            for c1 in range(channels):
                for c2 in range(channels):
                    covariance[k, c1, c2] = np.float32(covariance[k, c1, c2] / denominator)
    return covariance


@njit(cache=True, fastmath=False)
def _gaussian_memberships_kernel(values, mem, active, means, cov, priors, smooth):
    classes = mem.shape[3]
    determinants = np.abs(cov[:, 0, 0]).astype(np.float32)
    inverse = np.asarray(1.0 / cov[:, 0, 0], dtype=np.float32)
    distance = np.float32(0.0)
    total_energy = np.float32(0.0)
    pi_source = np.float32(3.14159265358979)
    log2pi = np.float32(np.log(np.float64(pi_source)))
    old = np.empty(classes, dtype=np.float32)
    for x in range(values.shape[0]):
        for y in range(values.shape[1]):
            for z in range(values.shape[2]):
                if active[x, y, z]:
                    for k in range(classes):
                        old[k] = mem[x, y, z, k]
                    denominator = np.float32(1e-30)
                    energy_terms = np.zeros(classes, dtype=np.float32)
                    for k in range(classes):
                        prior = priors[x, y, z, k]
                        if prior != 0.0:
                            normalizer = np.float32(-np.float32(0.5) * log2pi)
                            likelihood = (float(normalizer) - 0.5 * np.log(float(determinants[k])) + np.log(float(prior)))
                            energy_terms[k] = np.float32(likelihood)
                        if prior > 0.0:
                            delta = np.float32(values[x, y, z, 0] - means[0, k])
                            exponent = np.float32(np.float32(inverse[k] * delta) * delta)
                            exponent = np.float32(exponent / np.float32(2.0))
                            if smooth > 0.0:
                                neighbor_sum = np.float32(0.0)
                                for m in range(classes):
                                    if m != k:
                                        neighbor_sum = np.float32(neighbor_sum + mem[x+1, y, z, m])
                                for m in range(classes):
                                    if m != k:
                                        neighbor_sum = np.float32(neighbor_sum + mem[x-1, y, z, m])
                                for m in range(classes):
                                    if m != k:
                                        neighbor_sum = np.float32(neighbor_sum + mem[x, y+1, z, m])
                                for m in range(classes):
                                    if m != k:
                                        neighbor_sum = np.float32(neighbor_sum + mem[x, y-1, z, m])
                                for m in range(classes):
                                    if m != k:
                                        neighbor_sum = np.float32(neighbor_sum + mem[x, y, z+1, m])
                                for m in range(classes):
                                    if m != k:
                                        neighbor_sum = np.float32(neighbor_sum + mem[x, y, z-1, m])
                                exponent = np.float32(exponent + np.float32(smooth * neighbor_sum))
                            energy_terms[k] = np.float32(energy_terms[k] - exponent)
                            raw = (np.exp(np.float64(-exponent)) * np.float64(prior) / np.sqrt(np.float64(determinants[k])))
                            mem[x, y, z, k] = np.float32(raw)
                            denominator = np.float32(denominator + mem[x, y, z, k])
                        else:
                            mem[x, y, z, k] = np.float32(0.0)
                    if denominator > 0.0:
                        for k in range(classes):
                            mem[x, y, z, k] = np.float32(mem[x, y, z, k] / denominator)
                    else:
                        for k in range(classes):
                            mem[x, y, z, k] = np.float32(0.0)
                    for k in range(classes):
                        change = np.float32(abs(np.float32(mem[x, y, z, k] - old[k])))
                        if change > distance:
                            distance = change
                        contribution = float(energy_terms[k]) + np.log(float(mem[x, y, z, k])) * float(old[k])
                        total_energy = np.float32(float(total_energy) + contribution)
                else:
                    for k in range(classes):
                        mem[x, y, z, k] = np.float32(0.0)
    return mem, np.asarray([total_energy, distance], dtype=np.float32)


def compute_atlas_em_mixing_coefficients(
    memberships: ArrayLike, mask: ArrayLike, *, mode: str = "Class dependent"
) -> NDArray[np.float32]:
    """Reproduce source ``computePriors`` in x/y/z then class order."""
    mem = np.asarray(memberships, dtype=np.float32)
    active = np.asarray(mask, dtype=bool)
    if mem.ndim != 4 or any(n < 1 for n in mem.shape) or active.shape != mem.shape[:3]:
        raise ValueError("memberships must be (x, y, z, classes) with a matching nonempty mask")
    if not np.isfinite(mem).all():
        raise ValueError("memberships must be finite")
    if mode == "Equal":
        return np.full(mem.shape[-1], np.float32(1.0) / np.float32(mem.shape[-1]), dtype=np.float32)

    coefficients = np.zeros(mem.shape[-1], dtype=np.float32)
    denominator = np.float32(0.0)
    for x in range(mem.shape[0]):
        for y in range(mem.shape[1]):
            for z in range(mem.shape[2]):
                if active[x, y, z]:
                    denominator = np.float32(denominator + np.float32(1.0))
                    for k in range(mem.shape[-1]):
                        coefficients[k] = np.float32(coefficients[k] + mem[x, y, z, k])
    if denominator != np.float32(0.0):
        coefficients = np.asarray(coefficients / denominator, dtype=np.float32)
    return coefficients


def normalize_atlas_em_priors(atlas_priors: ArrayLike, mask: ArrayLike) -> NDArray[np.float32]:
    """Normalize transformed Atlas-mode class priors on masked voxels.

    This reproduces the normalization loop at the end of
    ``AtlasEMSegmentation.computeAtlas`` after atlas resampling. A zero-sum
    voxel is left as-is, and values outside the mask remain untouched.
    """
    priors = np.asarray(atlas_priors, dtype=np.float32)
    active = np.asarray(mask, dtype=bool)
    if priors.ndim != 4 or any(n < 1 for n in priors.shape[:3]) or priors.shape[3] < 1:
        raise ValueError("atlas_priors must have shape (x, y, z, classes)")
    if active.shape != priors.shape[:3]:
        raise ValueError("mask must match the atlas prior grid")
    if not np.isfinite(priors).all():
        raise ValueError("atlas priors must be finite")
    result = np.array(priors, dtype=np.float32, copy=True)
    for x in range(priors.shape[0]):
        for y in range(priors.shape[1]):
            for z in range(priors.shape[2]):
                if active[x, y, z]:
                    norm = np.float32(0.0)
                    for k in range(priors.shape[3]):
                        norm = np.float32(norm + result[x, y, z, k])
                    if norm != 0.0:
                        for k in range(priors.shape[3]):
                            result[x, y, z, k] = np.float32(result[x, y, z, k] / norm)
    return result


def compute_atlas_em_centroids(
    images: ArrayLike,
    memberships: ArrayLike,
    mask: ArrayLike,
    *,
    fields: ArrayLike | None = None,
    initial: ArrayLike | None = None,
) -> NDArray[np.float32]:
    """Reproduce source ``computeCentriods``; empty classes keep initial means."""
    values = np.asarray(images, dtype=np.float32)
    mem = np.asarray(memberships, dtype=np.float32)
    active = np.asarray(mask, dtype=bool)
    if values.ndim != 4 or any(n < 1 for n in values.shape):
        raise ValueError("images must have shape (x, y, z, channels)")
    if mem.ndim != 4 or mem.shape[:3] != values.shape[:3] or mem.shape[3] < 1:
        raise ValueError("memberships must match the image grid and contain classes")
    if active.shape != values.shape[:3]:
        raise ValueError("mask must match the image grid")
    field = None if fields is None else np.asarray(fields, dtype=np.float32)
    if field is not None and field.shape != values.shape:
        raise ValueError("fields must match image shape")
    if not np.isfinite(values).all() or not np.isfinite(mem).all() or (field is not None and not np.isfinite(field).all()):
        raise ValueError("images, memberships, and fields must be finite")
    means = np.zeros((values.shape[3], mem.shape[3]), dtype=np.float32) if initial is None else np.array(initial, dtype=np.float32, copy=True)
    if means.shape != (values.shape[3], mem.shape[3]):
        raise ValueError("initial centroids must have shape (channels, classes)")
    field_values = np.empty((0, 0, 0, 0), dtype=np.float32) if field is None else field
    return _centroids_kernel(values, mem, active, field_values, field is not None, means)


def compute_atlas_em_covariances(
    images: ArrayLike,
    memberships: ArrayLike,
    mask: ArrayLike,
    centroids: ArrayLike,
    *,
    fields: ArrayLike | None = None,
) -> NDArray[np.float32]:
    """Reproduce selected source covariance loop, including its denominator order.

    The Java source increments one shared class denominator inside both channel
    loops. That means it is incremented ``channels**2`` times per voxel; this
    source behavior is retained intentionally for direct source comparison.
    """
    values = np.asarray(images, dtype=np.float32)
    mem = np.asarray(memberships, dtype=np.float32)
    active = np.asarray(mask, dtype=bool)
    means = np.asarray(centroids, dtype=np.float32)
    if values.ndim != 4 or any(n < 1 for n in values.shape):
        raise ValueError("images must have shape (x, y, z, channels)")
    if mem.ndim != 4 or mem.shape[:3] != values.shape[:3] or mem.shape[3] < 1:
        raise ValueError("memberships must match the image grid and contain classes")
    channels, classes = values.shape[3], mem.shape[3]
    if active.shape != values.shape[:3] or means.shape != (channels, classes):
        raise ValueError("mask or centroids do not match the image grid and channels/classes")
    field = None if fields is None else np.asarray(fields, dtype=np.float32)
    if field is not None and field.shape != values.shape:
        raise ValueError("fields must match image shape")
    if not np.isfinite(values).all() or not np.isfinite(mem).all() or not np.isfinite(means).all() or (field is not None and not np.isfinite(field).all()):
        raise ValueError("images, memberships, centroids, and fields must be finite")
    field_values = np.empty((0, 0, 0, 0), dtype=np.float32) if field is None else field
    return _covariances_kernel(values, mem, active, means, field_values, field is not None)


def compute_atlas_em_atlas_prior_memberships(
    images: ArrayLike,
    memberships: ArrayLike,
    mask: ArrayLike,
    centroids: ArrayLike,
    covariance: ArrayLike,
    atlas_priors: ArrayLike,
    *,
    smoothing: float,
) -> tuple[NDArray[np.float32], NDArray[np.int8], np.float32]:
    """Apply one selected single-channel Atlas-prior EM membership pass.

    This is the wrapper's ``priorMode='Atlas'`` fixed-state branch, with source
    x/y/z and class iteration order and in-place six-neighbor smoothing. It
    returns updated memberships, selected hard class indices (-1 outside
    positive support), and maximum membership change. It does not deform the
    atlas or implement the full AtlasEM stopping/output wrapper.
    """
    values = np.asarray(images, dtype=np.float32)
    mem = np.array(memberships, dtype=np.float32, copy=True)
    active = np.asarray(mask, dtype=bool)
    means = np.asarray(centroids, dtype=np.float32)
    cov = np.asarray(covariance, dtype=np.float32)
    priors = np.asarray(atlas_priors, dtype=np.float32)
    if values.ndim != 4 or values.shape[3] != 1 or any(n < 3 for n in values.shape[:3]):
        raise ValueError("images must be a padded (x, y, z, 1) T1 grid with dimensions >= 3")
    if mem.ndim != 4 or mem.shape[:3] != values.shape[:3] or mem.shape[3] < 1:
        raise ValueError("memberships must match images and contain classes")
    classes = mem.shape[3]
    if active.shape != values.shape[:3] or priors.shape != mem.shape:
        raise ValueError("mask and atlas priors must match the image/membership grid")
    if means.shape != (1, classes) or cov.shape != (classes, 1, 1):
        raise ValueError("single-channel centroids/covariances must have one entry per class")
    if not np.isfinite(values).all() or not np.isfinite(mem).all() or not np.isfinite(means).all() or not np.isfinite(cov).all() or not np.isfinite(priors).all():
        raise ValueError("images, memberships, model parameters, and atlas priors must be finite")
    smooth = np.float32(smoothing)
    if not np.isfinite(smooth) or smooth < 0:
        raise ValueError("smoothing must be finite and nonnegative")
    determinants = np.abs(cov[:, 0, 0]).astype(np.float32)
    if np.any(determinants == 0):
        raise ValueError("selected AtlasEM covariance must be nonsingular")
    inverses = np.reciprocal(cov[:, 0, 0], dtype=np.float32)
    distance = np.float32(0.0)

    for x in range(values.shape[0]):
        for y in range(values.shape[1]):
            for z in range(values.shape[2]):
                if not active[x, y, z]:
                    mem[x, y, z, :] = 0.0
                    continue
                old = mem[x, y, z, :].copy()
                denominator = np.float32(1.0e-30)
                for k in range(classes):
                    prior = priors[x, y, z, k]
                    if prior > 0.0:
                        difference = np.float32(values[x, y, z, 0] - means[0, k])
                        exponent = np.float32(inverses[k] * difference)
                        exponent = np.float32(exponent * difference)
                        exponent = np.float32(exponent / np.float32(2.0))
                        if smooth > 0.0:
                            neighbor_sum = np.float32(0.0)
                            for nx, ny, nz in ((x + 1, y, z), (x - 1, y, z), (x, y + 1, z),
                                               (x, y - 1, z), (x, y, z + 1), (x, y, z - 1)):
                                for m in range(classes):
                                    if m != k:
                                        neighbor_sum = np.float32(neighbor_sum + mem[nx, ny, nz, m])
                            exponent = np.float32(exponent + np.float32(smooth * neighbor_sum))
                        raw = np.exp(np.float64(-exponent)) * np.float64(prior) / np.sqrt(np.float64(determinants[k]))
                        mem[x, y, z, k] = np.float32(raw)
                        denominator = np.float32(denominator + mem[x, y, z, k])
                    else:
                        mem[x, y, z, k] = 0.0
                if denominator > 0.0:
                    for k in range(classes):
                        mem[x, y, z, k] = np.float32(mem[x, y, z, k] / denominator)
                else:
                    mem[x, y, z, :] = 0.0
                for k in range(classes):
                    change = np.float32(abs(np.float32(mem[x, y, z, k] - old[k])))
                    if change > distance:
                        distance = change

    labels = np.full(values.shape[:3], -1, dtype=np.int8)
    for x in range(values.shape[0]):
        for y in range(values.shape[1]):
            for z in range(values.shape[2]):
                best = 0
                for k in range(1, classes):
                    if mem[x, y, z, k] > mem[x, y, z, best]:
                        best = k
                if mem[x, y, z, best] > 0:
                    labels[x, y, z] = best
    return mem, labels, distance


def compute_atlas_em_diagonal_covariances(
    images: ArrayLike,
    memberships: ArrayLike,
    centroids: ArrayLike,
    *,
    fields: ArrayLike | None = None,
) -> NDArray[np.float32]:
    """Reproduce ``computeDiagonalCovariances``' standard-deviation estimates.

    The source adds ``membership² + 1e-30f`` to numerator and denominator for
    every voxel, including outside the constructor's mask. The returned
    standard deviations are installed directly as diagonal covariance values
    by the selected wrapper, preserving that historical behavior.
    """
    values = np.asarray(images, dtype=np.float32)
    mem = np.asarray(memberships, dtype=np.float32)
    if values.ndim != 4 or any(n < 1 for n in values.shape) or values.shape[3] < 1:
        raise ValueError("images must have shape (x, y, z, channels)")
    if mem.ndim != 4 or mem.shape[:3] != values.shape[:3] or mem.shape[3] < 1:
        raise ValueError("memberships must match the image grid and contain classes")
    channels, classes = values.shape[3], mem.shape[3]
    means = np.asarray(centroids, dtype=np.float32)
    if means.shape != (channels, classes):
        raise ValueError("centroids must have shape (channels, classes)")
    field = None if fields is None else np.asarray(fields, dtype=np.float32)
    if field is not None and field.shape != values.shape:
        raise ValueError("fields must match image shape")
    if not np.isfinite(values).all() or not np.isfinite(mem).all() or not np.isfinite(means).all() or (field is not None and not np.isfinite(field).all()):
        raise ValueError("images, memberships, centroids, and fields must be finite")

    result = np.zeros((channels, classes), dtype=np.float32)
    epsilon = np.float32(1.0e-30)
    for c in range(channels):
        for k in range(classes):
            numerator = np.float32(0.0)
            denominator = np.float32(0.0)
            for x in range(values.shape[0]):
                for y in range(values.shape[1]):
                    for z in range(values.shape[2]):
                        membership = mem[x, y, z, k]
                        weight = np.float32(np.float32(membership * membership) + epsilon)
                        if field is None:
                            difference = np.float32(values[x, y, z, c] - means[c, k])
                        else:
                            corrected = np.float32(field[x, y, z, c] * values[x, y, z, c])
                            difference = np.float32(corrected - means[c, k])
                        term = np.float32(weight * difference)
                        term = np.float32(term * difference)
                        numerator = np.float32(numerator + term)
                        denominator = np.float32(denominator + weight)
            result[c, k] = np.float32(np.sqrt(np.float64(numerator / denominator))) if denominator > 0.0 else np.float32(0.0)
    return result


def update_atlas_em_gaussian_memberships(
    images: ArrayLike,
    memberships: ArrayLike,
    mask: ArrayLike,
    centroids: ArrayLike,
    covariance: ArrayLike,
    atlas_priors: ArrayLike,
    *,
    smoothing: float = 0.0,
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Perform one source-order Atlas-mode Gaussian membership update.

    This is the selected one-channel/no-field branch of
    ``AtlasEMSegmentation.computeMemberships``. Updates occur in place in
    x/y/z/class order, so later voxels' six-neighbor smoothing sees already
    updated memberships. The returned pair is ``(memberships, [energy,
    maximum_change])``. As in the Java method, energy arithmetic can produce
    NaN when a zero old membership is multiplied by ``log(0)``; callers that
    compare the source stop metric must preserve and report that behavior.
    The padded grid must have a false mask on all six boundary faces.
    """
    values = np.asarray(images, dtype=np.float32)
    mem = np.array(memberships, dtype=np.float32, copy=True)
    active = np.asarray(mask, dtype=bool)
    means = np.asarray(centroids, dtype=np.float32)
    cov = np.asarray(covariance, dtype=np.float32)
    priors = np.asarray(atlas_priors, dtype=np.float32)
    if values.ndim != 4 or values.shape[3] != 1 or any(n < 3 for n in values.shape[:3]):
        raise ValueError("images must be a padded (x, y, z, 1) single-channel array")
    shape = values.shape[:3]
    classes = mem.shape[-1] if mem.ndim == 4 else 0
    if mem.ndim != 4 or mem.shape[:3] != shape or classes < 1:
        raise ValueError("memberships must match the image grid and contain classes")
    if active.shape != shape or means.shape != (1, classes) or cov.shape != (classes, 1, 1) or priors.shape != mem.shape:
        raise ValueError("mask, centroids, covariance, or atlas priors have incompatible shapes")
    if not all(np.isfinite(a).all() for a in (values, mem, means, cov, priors)):
        raise ValueError("images, memberships, centroids, covariance, and priors must be finite")
    if np.any(priors < 0) or np.any(cov <= 0) or not np.isfinite(smoothing) or smoothing < 0:
        bad_covariances = np.flatnonzero(cov[:, 0, 0] <= 0).tolist()
        negative_prior_count = int(np.count_nonzero(priors < 0))
        raise ValueError(
            "Atlas priors and smoothing must be nonnegative; covariance must be positive "
            f"(nonpositive covariance classes={bad_covariances}, "
            f"negative prior voxels={negative_prior_count}, "
            f"covariance diagonal={cov[:, 0, 0].tolist()})"
        )
    if (active[0].any() or active[-1].any() or active[:, 0].any() or active[:, -1].any()
            or active[:, :, 0].any() or active[:, :, -1].any()):
        raise ValueError("mask must be false on all boundary faces to provide the source padding")

    smooth = np.float32(smoothing)
    return _gaussian_memberships_kernel(values, mem, active, means, cov, priors, smooth)
