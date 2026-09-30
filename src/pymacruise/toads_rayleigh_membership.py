"""One fixed-input homogeneous RAYLEIGH membership pass for base TOADS.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``. This
ports a selected class-energy/membership slice with fixed parameters; it does
not establish the complete atlas-driven classification procedure.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_math import toads_inverse_shape_prior_energy, toads_ordered_membership_sweep
from .toads_rayleigh_energy import toads_rayleigh_intensity_energy


@dataclass(frozen=True)
class ToadsRayleighMembershipPass:
    """Energy map, updated memberships, and objective from one sweep."""

    data_energy: NDArray[np.float32]
    memberships: NDArray[np.float32]
    objective: np.float32


def toads_rayleigh_membership_pass(
    image_values: ArrayLike,
    initial_memberships: ArrayLike,
    centroids: ArrayLike,
    modality_weights: ArrayLike,
    class_variance_weights: ArrayLike,
    image_variances: ArrayLike,
    transformed_shapes: ArrayLike,
    closest_labels: ArrayLike,
    smoothing_priors: ArrayLike,
    centroid_similarity: ArrayLike,
    atlas_shape_priors: ArrayLike,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    *,
    intensity_field: ArrayLike | None = None,
    psf_bias: ArrayLike | None = None,
) -> ToadsRayleighMembershipPass:
    """Compose RAYLEIGH data energy, atlas prior, and the shared TOADS sweep.

    Every class is assumed to be RAYLEIGH. Centroids are the class/channel
    Rayleigh scales and must be finite and positive so the shared sweep
    receives finite energies. All image statistics, shape maps, relations,
    similarities, and prior weights are explicit fixed inputs; this function
    does not estimate them or implement TOADS/MGDM topology and outer loops.
    Arrays are already aligned in true-RAS voxel order.
    """

    images = np.asarray(image_values, dtype=np.float32)
    initial = np.asarray(initial_memberships, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    if images.ndim != 4 or initial.ndim != 4 or images.shape[:3] != initial.shape[:3]:
        raise ValueError("image_values and initial_memberships must be 4-D on the same grid")
    if centers.shape != (images.shape[-1], initial.shape[-1]):
        raise ValueError("centroids must have shape (modality,class)")
    if not np.isfinite(centers).all() or np.any(centers <= np.float32(0.0)):
        raise ValueError("RAYLEIGH scales must be finite and positive")

    data_energy = toads_rayleigh_intensity_energy(
        images,
        centers,
        modality_weights,
        class_variance_weights,
        image_variances,
        intensity_field=intensity_field,
        psf_bias=psf_bias,
    )
    inverse_prior, prior_objective = toads_inverse_shape_prior_energy(
        transformed_shapes,
        atlas_shape_priors,
        centroid_similarity,
        prior_coefficient,
    )
    updated, objective = toads_ordered_membership_sweep(
        data_energy,
        initial,
        transformed_shapes,
        closest_labels,
        smoothing_priors,
        smoothing,
        prior_factor,
        prior_inverse_energy=inverse_prior,
        prior_objective_energy=prior_objective,
    )
    return ToadsRayleighMembershipPass(data_energy, updated, objective)
