"""Source-ordered fixed-prior AtlasEM iterations for synthetic/reference use."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.atlas_em_statistics import (
    compute_atlas_em_centroids,
    compute_atlas_em_covariances,
    update_atlas_em_gaussian_memberships,
)


@dataclass(frozen=True)
class AtlasEMIterationResult:
    memberships: NDArray[np.float32]
    centroids: NDArray[np.float32]
    covariances: NDArray[np.float32]
    energy_history: NDArray[np.float32]
    maximum_change_history: NDArray[np.float32]
    iterations: int
    converged: bool
    regularized_classes: tuple[int, ...] = ()


def run_atlas_em_fixed_prior_iterations(
    images: ArrayLike,
    memberships: ArrayLike,
    mask: ArrayLike,
    centroids: ArrayLike,
    covariances: ArrayLike,
    atlas_priors: ArrayLike,
    *,
    max_iterations: int = 99,
    max_membership_change: float = 0.01,
    smoothing: float = 10.0,
    covariance_floor: float = 1.0e-6,
) -> AtlasEMIterationResult:
    """Run the selected no-field, one-channel Atlas-mode outer EM loop.

    This assembles the source wrapper loop after AtlasEM initialization and
    transformed-prior generation: update centroids, update covariance, update
    posteriors, then stop on the source maximum-membership-change threshold or
    iteration cap. Priors are fixed across iterations, matching the selected
    `AlgorithmAtlasEMSegmentation` loop where `computeAtlas()` is not called
    again. It deliberately accepts the (potentially synthetic) transformed
    prior as input and does not perform registration or Dura evolution.
    """
    values = np.asarray(images, dtype=np.float32)
    mem = np.array(memberships, dtype=np.float32, copy=True)
    active = np.asarray(mask, dtype=bool)
    means = np.asarray(centroids, dtype=np.float32).copy()
    cov = np.asarray(covariances, dtype=np.float32).copy()
    priors = np.asarray(atlas_priors, dtype=np.float32)
    if values.ndim != 4 or values.shape[3] != 1 or any(n < 3 for n in values.shape[:3]):
        raise ValueError("images must be a padded (x, y, z, 1) grid")
    if mem.ndim != 4 or mem.shape[:3] != values.shape[:3] or mem.shape[3] < 1:
        raise ValueError("memberships must match images and contain classes")
    classes = mem.shape[3]
    if active.shape != values.shape[:3] or priors.shape != mem.shape:
        raise ValueError("mask and atlas priors must match the image/membership grid")
    if means.shape != (1, classes) or cov.shape != (classes, 1, 1):
        raise ValueError("centroids/covariances must match single-channel classes")
    if not all(np.isfinite(a).all() for a in (values, mem, means, cov, priors)):
        raise ValueError("images, memberships, model parameters, and priors must be finite")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 0:
        raise ValueError("max_iterations must be a nonnegative integer")
    if not np.isfinite(max_membership_change) or max_membership_change < 0:
        raise ValueError("max_membership_change must be finite and nonnegative")
    if not np.isfinite(covariance_floor) or covariance_floor <= 0:
        raise ValueError("covariance_floor must be finite and positive")

    energies: list[np.float32] = []
    changes: list[np.float32] = []
    regularized: set[int] = set()
    converged = False
    for _ in range(max_iterations):
        means = compute_atlas_em_centroids(values, mem, active)
        cov = compute_atlas_em_covariances(values, mem, active, means)
        floor = np.float32(covariance_floor)
        for class_index in range(classes):
            if cov[class_index, 0, 0] < floor:
                cov[class_index, 0, 0] = floor
                regularized.add(class_index)
        mem, metrics = update_atlas_em_gaussian_memberships(
            values, mem, active, means, cov, priors, smoothing=smoothing
        )
        energies.append(metrics[0])
        changes.append(metrics[1])
        if metrics[1] < np.float32(max_membership_change):
            converged = True
            break

    return AtlasEMIterationResult(
        memberships=mem,
        centroids=means,
        covariances=cov,
        energy_history=np.asarray(energies, dtype=np.float32),
        maximum_change_history=np.asarray(changes, dtype=np.float32),
        iterations=len(changes),
        converged=converged,
        regularized_classes=tuple(sorted(regularized)),
    )
