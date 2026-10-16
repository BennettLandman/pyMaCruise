"""Selected fixed-input AtlasEM statistics updates.

This module reproduces the active ``AtlasEMSegmentation.computePriors``,
``computeCentriods`` (source spelling), and ``computeCovariance`` loops. It is
not the AtlasEM pipeline: atlas registration, preprocessing, membership updates,
and the missing selected atlas files are tracked under P44b.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


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
                            if field is not None:
                                weighted = np.float32(membership * field[x, y, z, c])
                                weighted = np.float32(weighted * sample)
                            else:
                                weighted = np.float32(membership * sample)
                            numerator = np.float32(numerator + weighted)
                            denominator = np.float32(denominator + mem[x, y, z, k])
            if denominator != np.float32(0.0):
                means[c, k] = np.float32(numerator / denominator)
    return means


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
    covariance = np.zeros((classes, channels, channels), dtype=np.float32)
    for k in range(classes):
        denominator = np.float32(0.0)
        for x in range(values.shape[0]):
            for y in range(values.shape[1]):
                for z in range(values.shape[2]):
                    if not active[x, y, z]:
                        continue
                    membership = mem[x, y, z, k]
                    samples = values[x, y, z]
                    if field is not None:
                        samples = np.asarray(field[x, y, z] * samples, dtype=np.float32)
                    difference = np.asarray(samples - means[:, k], dtype=np.float32)
                    for c1 in range(channels):
                        for c2 in range(channels):
                            product = np.float32(difference[c1] * difference[c2])
                            product = np.float32(product * membership)
                            covariance[k, c1, c2] = np.float32(covariance[k, c1, c2] + product)
                            denominator = np.float32(denominator + membership)
        if denominator != np.float32(0.0):
            for c1 in range(channels):
                for c2 in range(channels):
                    covariance[k, c1, c2] = np.float32(covariance[k, c1, c2] / denominator)
    return covariance


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
