"""Compose the selected AtlasEM registration result with its fixed-prior EM tail."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from pymacruise.atlas_em_iteration import (
    AtlasEMIterationResult,
    run_atlas_em_fixed_prior_iterations,
)
from pymacruise.atlas_em_preprocessing import (
    build_atlas_em_internal_mask,
    transform_atlas_em_rigid_priors,
)
from pymacruise.atlas_em_registration import AtlasEMTwoStageRegistrationResult
from pymacruise.atlas_em_statistics import (
    compute_atlas_em_atlas_prior_memberships,
    compute_atlas_em_diagonal_covariances,
)
from pymacruise.atlas_em_outputs import (
    generate_atlas_em_dura_inputs,
    generate_atlas_em_segmentation,
)
from pymacruise.atlas_em_descriptor import parse_atlas_em_descriptor


@dataclass(frozen=True)
class AtlasEMRegisteredExecutionResult:
    """Transformed priors, EM results, tissue indices and Dura inputs."""

    mask: NDArray[np.bool_]
    priors: NDArray[np.float32]
    initial_memberships: NDArray[np.float32]
    initial_metrics: NDArray[np.float32]
    iterations: AtlasEMIterationResult
    segmentation: NDArray[np.int8]
    dura_inputs: NDArray[np.float32]


def finish_atlas_em_registered_execution(
    registration: AtlasEMTwoStageRegistrationResult,
    descriptor_path: str,
    *,
    max_iterations: int = 99,
    max_membership_change: float = 0.01,
    smoothing: float = 10.0,
) -> AtlasEMRegisteredExecutionResult:
    """Run AtlasEM's initialization membership update, fixed-prior loop, outputs.

    The selected wrapper computes diagonal covariance estimates, transforms
    and normalizes the registered shape priors, performs an initial Atlas-mode
    membership update, then iterates centroid/covariance/posterior updates
    until the source maximum-change threshold or cap. This function consumes
    the preceding two-stage descriptor registration output. It does not run
    the pre-LesionTOADS Dura level-set evolution or export NIfTI files.
    """
    images = np.asarray(registration.recropped.images, dtype=np.float32)
    memberships = np.asarray(registration.post_recrop_initialization.memberships, dtype=np.float32)
    centroids = np.asarray(registration.post_recrop_initialization.centroids, dtype=np.float32)
    if images.ndim != 4 or images.shape[-1] != 1:
        raise ValueError("selected AtlasEM execution currently requires one T1 channel")
    if memberships.shape[:3] != images.shape[:3] or centroids.shape != (1, memberships.shape[-1]):
        raise ValueError("registered AtlasEM initialization does not match the recropped T1")
    descriptor = parse_atlas_em_descriptor(descriptor_path)
    if len(descriptor.structures) != memberships.shape[-1]:
        raise ValueError("descriptor class count does not match registered memberships")

    active = build_atlas_em_internal_mask(images)
    dimensions = tuple(int(value) for value in images.shape[:3])
    priors = transform_atlas_em_rigid_priors(
        registration.atlas.shape_priors,
        active,
        image_dimensions_xyz=dimensions,
        image_resolutions_mm=registration.recropped.signed_resolutions_mm,
        template_resolutions_mm=descriptor.shape_resolutions_mm,
        transform=registration.post_recrop_registration.transform,
        shape_slope=descriptor.shape_slope,
    )

    initial_deviations = compute_atlas_em_diagonal_covariances(
        images, memberships, centroids
    )
    covariance = np.zeros((memberships.shape[-1], 1, 1), dtype=np.float32)
    covariance[:, 0, 0] = initial_deviations[0]
    initialized_memberships, _initial_labels, initial_metrics = compute_atlas_em_atlas_prior_memberships(
        images,
        memberships,
        active,
        centroids,
        covariance,
        priors,
        smoothing=smoothing,
    )
    iterations = run_atlas_em_fixed_prior_iterations(
        images,
        initialized_memberships,
        active,
        centroids,
        covariance,
        priors,
        max_iterations=max_iterations,
        max_membership_change=max_membership_change,
        smoothing=smoothing,
    )
    segmentation = generate_atlas_em_segmentation(iterations.memberships)
    dura_inputs = generate_atlas_em_dura_inputs(
        iterations.memberships,
        segmentation,
        active,
        tuple(item.name for item in descriptor.structures),
    )
    return AtlasEMRegisteredExecutionResult(
        active,
        priors,
        initialized_memberships,
        np.asarray(initial_metrics, dtype=np.float32),
        iterations,
        segmentation,
        dura_inputs,
    )
