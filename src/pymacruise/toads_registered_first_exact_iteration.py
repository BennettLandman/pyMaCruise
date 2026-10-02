"""Continue registered unity-field startup through the first loop membership."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_default_atlas_membership import (
    ToadsDefaultAtlasMembershipPass,
    toads_default_atlas_membership_pass,
)
from .toads_default_exact_centroids import toads_default_exact_centroid_update
from .toads_default_exact_variances import toads_default_exact_class_variances
from .toads_inhomogeneity_transition import toads_next_inhomogeneity_degree
from .toads_math import toads_optimized_centroid_similarity, toads_optimized_image_variance
from .toads_registered_unity_field_startup import ToadsRegisteredUnityFieldStartup


@dataclass(frozen=True)
class ToadsRegisteredFirstExactIteration:
    """Outputs of the first exact update after startup thinning."""

    startup: ToadsRegisteredUnityFieldStartup
    degree: int
    fields: NDArray[np.float32]
    centroids: NDArray[np.float32]
    class_standard_deviations: NDArray[np.float32]
    image_variance_weights: NDArray[np.float32]
    image_standard_deviation: NDArray[np.float32]
    centroid_similarity: NDArray[np.float32]
    membership_pass: ToadsDefaultAtlasMembershipPass
    stage_order: tuple[str, ...]


def continue_registered_unity_startup_through_first_exact_iteration(
    alignment_run: Any,
    startup: ToadsRegisteredUnityFieldStartup,
    *,
    category_exponents: Any,
    exponent_by_channel_class: Any,
    modality_weights: Any,
    class_variance_weights: Any,
    smoothing_priors: Any,
    atlas_shape_priors: Any,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    prior_scale: float = 0.05,
    low_energy_ratio: float = 0.01,
    polynomial_degree: int = 3,
    centroid_smoothness: float = 0.1,
    psf_bias: Any | None = None,
) -> ToadsRegisteredFirstExactIteration:
    """Run the first loop's degree-zero field no-op and exact update stages.

    The input is the selected registered startup after its initial exact
    membership and first thinning. Initial ``Eratio=1`` and degree zero do
    not schedule a field fit under the selected default threshold. Then the
    wrapper refreshes centroids/deviations and executes the next exact
    membership pass. Transform, shape and closest-label maps come from the
    supplied alignment run and remain fixed at this boundary.
    """
    if not isinstance(startup, ToadsRegisteredUnityFieldStartup):
        raise TypeError("startup must be a registered unity-field startup result")
    degree = toads_next_inhomogeneity_degree(
        startup.degree_state,
        polynomial_degree,
        startup.energy_ratio,
        low_energy_ratio,
    )
    if degree != 0:
        raise ValueError("the selected initial energy ratio must keep the first loop field degree at zero")
    try:
        images = np.moveaxis(
            np.asarray(alignment_run.template_refresh.recropped_images, dtype=np.float32),
            0,
            -1,
        )
        previous_centroids = np.asarray(
            alignment_run.segmentation_initialization.centroids, dtype=np.float32
        )
        closest = np.asarray(alignment_run.relation_state.closest_labels)
        shapes = np.asarray(alignment_run.transformed_shape_priors, dtype=np.float32)
    except AttributeError as exc:
        raise TypeError("alignment_run is missing selected registered TOADS state") from exc

    memberships = np.asarray(
        startup.statistics_membership.membership_pass.memberships, dtype=np.float32
    )
    labels = np.asarray(startup.thinning.segmentation)
    fields = startup.fields
    if images.ndim != 4 or images.shape != fields.shape:
        raise ValueError("registered images and unity fields must share channel-last dimensions")
    channels = images.shape[-1]
    if previous_centroids.shape != (channels, 11):
        raise ValueError("registered centroids must match image channels and default atlas classes")
    if memberships.shape != images.shape[:3] + (11,) or labels.shape != images.shape[:3]:
        raise ValueError("startup memberships and thinned labels must match the registered image grid")
    if shapes.shape != memberships.shape or closest.shape != images.shape[:3] + (4,):
        raise ValueError("registered shape and closest-label maps must match startup state")

    centers = toads_default_exact_centroid_update(
        images,
        memberships,
        previous_centroids,
        category_exponents,
        intensity_field=fields,
        smoothness=centroid_smoothness,
    )
    deviations = toads_default_exact_class_variances(
        images,
        memberships,
        centers,
        labels,
        intensity_field=fields,
    )
    image_variance, image_std = toads_optimized_image_variance(
        images,
        memberships[..., 1:],
        centers[:, 1:],
        intensity_field=fields,
    )
    similarity = toads_optimized_centroid_similarity(
        centers,
        deviations,
        image_variance,
        prior_scale,
    )
    membership = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=memberships,
        centroids=centers,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        image_variances=image_variance,
        transformed_shapes=shapes,
        closest_labels=closest,
        smoothing_priors=smoothing_priors,
        centroid_similarity=similarity,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        intensity_field=fields,
        psf_bias=psf_bias,
    )
    stages = (
        *startup.stage_order,
        "iteration1.computeCorrectionField.degree0.noop",
        "iteration1.computeExactCentroids",
        "iteration1.computeExactVariances",
        "iteration1.computeImageVariance",
        "iteration1.computeCentroidSimilarity",
        "iteration1.computeExactMemberships",
    )
    return ToadsRegisteredFirstExactIteration(
        startup=startup,
        degree=degree,
        fields=fields,
        centroids=centers,
        class_standard_deviations=deviations,
        image_variance_weights=image_variance,
        image_standard_deviation=image_std,
        centroid_similarity=similarity,
        membership_pass=membership,
        stage_order=stages,
    )
