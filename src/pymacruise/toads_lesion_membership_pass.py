"""Fixed-input base-TOADS OBJECT+LESION membership pass.

This composes the selected ``MultipleToadSegmentation.computeExactMemberships``
OBJECT and LESION branches. The source's categorical segmentation gate controls
both its second LESION-centroid energy and the separate lesion-membership map.
It is a limited selected-source pass, not a complete TOADS iteration.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_math import (
    toads_inverse_shape_prior_energy,
    toads_object_intensity_energy,
    toads_shape_relation_factors,
)


@dataclass(frozen=True)
class ToadsLesionMembershipPass:
    """Raw terms, objective energies, and updated OBJECT+LESION maps."""

    class_data_terms: NDArray[np.float32]
    lesion_data_terms: NDArray[np.float32]
    class_energy: NDArray[np.float32]
    memberships: NDArray[np.float32]
    lesion_membership: NDArray[np.float32]
    denominator: NDArray[np.float32]
    objective: np.float32


def toads_lesion_membership_pass(
    image_values: ArrayLike,
    initial_memberships: ArrayLike,
    initial_lesion_membership: ArrayLike,
    segmentation: ArrayLike,
    class_types: Sequence[str],
    centroids: ArrayLike,
    lesion_centroid: ArrayLike,
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
    zero_threshold: float = 1.0e-10,
    infinite_reciprocal: float = 1.0e10,
) -> ToadsLesionMembershipPass:
    """Run one source-ordered OBJECT+LESION pass with fixed inputs.

    Input arrays use RAS ``(x,y,z[,modality/class])`` indexing. The class axis
    must contain exactly one OBJECT and one LESION in either order. The two
    modality-weight columns are source rows 0 (ordinary residual) and 1
    (segmentation-gated lesion residual). Field multiplication occurs before
    PSF-bias and centroid subtraction.

    The pass updates interior voxels in x/y/z order and classes in their given
    order, in place. For the LESION class, raw ``numL`` exists only where the
    segmentation label equals that class. Shape-skipped classes and
    DISTANCE relations are not part of this exact-membership path.
    """
    images = np.asarray(image_values, dtype=np.float32)
    initial = np.asarray(initial_memberships, dtype=np.float32)
    lesion_initial = np.asarray(initial_lesion_membership, dtype=np.float32)
    labels = np.asarray(segmentation)
    types = tuple(str(value).upper() for value in class_types)
    centers = np.asarray(centroids, dtype=np.float32)
    lesion_center = np.asarray(lesion_centroid, dtype=np.float32)
    weights = np.asarray(modality_weights, dtype=np.float32)
    variance_weights = np.asarray(class_variance_weights, dtype=np.float32)
    variances = np.asarray(image_variances, dtype=np.float32)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    closest = np.asarray(closest_labels)
    smooth_priors = np.asarray(smoothing_priors, dtype=np.float32)
    similarities = np.asarray(centroid_similarity, dtype=np.float32)
    shape_priors = np.asarray(atlas_shape_priors, dtype=np.float32)

    if images.ndim != 4 or images.shape[-1] < 1 or min(images.shape[:3]) < 3:
        raise ValueError("image_values must have shape (x>=3,y>=3,z>=3,modalities>=1)")
    spatial = images.shape[:3]
    channels = images.shape[-1]
    classes = len(types)
    if types.count("OBJECT") != 1 or types.count("LESION") != 1 or classes != 2:
        raise ValueError("class_types must contain exactly one OBJECT and one LESION")
    if initial.shape != spatial + (classes,) or not np.isfinite(initial).all():
        raise ValueError("initial_memberships must be finite and match the two-class grid")
    if lesion_initial.shape != spatial or not np.isfinite(lesion_initial).all():
        raise ValueError("initial_lesion_membership must be finite and match the spatial grid")
    if labels.shape != spatial or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer map matching the spatial grid")
    if centers.shape != (channels, classes) or not np.isfinite(centers).all():
        raise ValueError("centroids must be finite with shape (modality, class)")
    if lesion_center.shape != (channels,) or not np.isfinite(lesion_center).all():
        raise ValueError("lesion_centroid must be finite with one value per modality")
    if weights.shape != (channels, 2) or not np.isfinite(weights).all():
        raise ValueError("modality_weights must provide source rows 0 and 1 per modality")
    if variance_weights.shape != (channels, classes) or not np.isfinite(variance_weights).all():
        raise ValueError("class_variance_weights must match (modality, class)")
    if variances.shape != (channels,) or not np.isfinite(variances).all():
        raise ValueError("image_variances must contain one finite value per modality")
    if shapes.shape != initial.shape or not np.isfinite(shapes).all():
        raise ValueError("transformed_shapes must be finite and match membership dimensions")
    if closest.shape != spatial + (4,) or not np.issubdtype(closest.dtype, np.integer):
        raise ValueError("closest_labels must be an integer four-slot map on the spatial grid")
    if smooth_priors.shape != (classes,) or not np.isfinite(smooth_priors).all():
        raise ValueError("smoothing_priors must contain one finite value per class")
    if similarities.shape != (classes, classes) or not np.isfinite(similarities).all():
        raise ValueError("centroid_similarity must be finite and square by class")
    if shape_priors.shape != (classes,) or not np.isfinite(shape_priors).all():
        raise ValueError("atlas_shape_priors must contain one finite value per class")
    if not np.isfinite(images).all():
        raise ValueError("image_values must be finite")

    field = None if intensity_field is None else np.asarray(intensity_field, dtype=np.float32)
    bias = None if psf_bias is None else np.asarray(psf_bias, dtype=np.float32)
    if field is not None and (field.shape != images.shape or not np.isfinite(field).all()):
        raise ValueError("intensity_field must be finite and match image_values")
    if bias is not None and (bias.shape != images.shape or not np.isfinite(bias).all()):
        raise ValueError("psf_bias must be finite and match image_values")

    smooth = np.float32(smoothing)
    prior_factor_value = np.float32(prior_factor)
    coefficient = np.float32(prior_coefficient)
    zero = np.float32(zero_threshold)
    infinity = np.float32(infinite_reciprocal)
    if not np.isfinite([smooth, prior_factor_value, coefficient, zero, infinity]).all():
        raise ValueError("smoothing, prior, and reciprocal parameters must be finite")
    if smooth < 0 or zero < 0 or infinity <= 0:
        raise ValueError("smoothing/threshold must be nonnegative and reciprocal cap positive")

    # `shape` is only populated in the source when priorCoefficient > 0.
    # With zero coefficient it remains the zero-initialized array; DIRECT-only
    # references can still validate smoothing independently of that quirk.
    relation_shapes = shapes if coefficient > np.float32(0.0) else np.zeros_like(shapes)
    factors, suppressed = toads_shape_relation_factors(
        relation_shapes, closest, float(prior_factor_value)
    )
    if np.any(suppressed[1:-1, 1:-1, 1:-1]):
        raise ValueError("all interior classes must have active shapes in this fixed-input slice")

    object_terms = toads_object_intensity_energy(
        images,
        centers,
        weights[:, 0],
        variance_weights,
        variances,
        intensity_field=field,
        psf_bias=bias,
    )
    lesion_index = types.index("LESION")
    lesion_terms = np.zeros(initial.shape, dtype=np.float32)
    corrected = images if field is None else np.asarray(field * images, dtype=np.float32)
    if bias is not None:
        corrected = np.asarray(corrected - bias, dtype=np.float32)
    for x in range(spatial[0]):
        for y in range(spatial[1]):
            for z in range(spatial[2]):
                if labels[x, y, z] != lesion_index:
                    continue
                value = np.float32(0.0)
                for channel in range(channels):
                    residual = np.float32(corrected[x, y, z, channel] - lesion_center[channel])
                    weight = np.float32(weights[channel, 1] * variance_weights[channel, lesion_index])
                    weight = np.float32(weight * variances[channel])
                    contribution = np.float32(weight * residual)
                    contribution = np.float32(contribution * residual)
                    value = np.float32(value + contribution)
                lesion_terms[x, y, z, lesion_index] = value

    prior_inverse, prior_objective = toads_inverse_shape_prior_energy(
        shapes, shape_priors, similarities, float(coefficient)
    )
    current = initial.copy()
    lesion = lesion_initial.copy()
    denominator_map = np.zeros(spatial, dtype=np.float32)
    energy_map = np.zeros(initial.shape, dtype=np.float32)
    objective = np.float32(0.0)
    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))

    for x in range(1, spatial[0] - 1):
        for y in range(1, spatial[1] - 1):
            for z in range(1, spatial[2] - 1):
                denominator = np.float32(0.0)
                objective_at_voxel = np.zeros(classes, dtype=np.float32)
                for k, kind in enumerate(types):
                    is_assigned_lesion = kind == "LESION" and labels[x, y, z] == k
                    num = np.float32(object_terms[x, y, z, k])
                    num_lesion = np.float32(lesion_terms[x, y, z, k]) if is_assigned_lesion else np.float32(0.0)
                    if is_assigned_lesion:
                        energy = np.float32(
                            np.float32(np.float32(1.0) - lesion[x, y, z]) * num
                            + np.float32(lesion[x, y, z] * num_lesion)
                        )
                    else:
                        energy = num

                    if smooth > np.float32(0.0):
                        neighbor_sum = np.float32(0.0)
                        for dx, dy, dz in offsets:
                            for other in range(classes):
                                if other == k:
                                    continue
                                value = current[x + dx, y + dy, z + dz, other]
                                neighbor_sum = np.float32(
                                    neighbor_sum + np.float32(value * value)
                                )
                        smooth_increment = np.float32(
                            np.float32(smooth_priors[k] * smooth) * neighbor_sum
                        )
                        num = np.float32(num + smooth_increment)
                        if is_assigned_lesion:
                            num_lesion = np.float32(num_lesion + smooth_increment)
                        smooth_objective = np.float32(
                            np.float32(np.float32(0.5) * smooth_priors[k]) * smooth
                        )
                        smooth_objective = np.float32(smooth_objective * neighbor_sum)
                        energy = np.float32(energy + smooth_objective)

                    if coefficient > np.float32(0.0):
                        prior_increment = prior_inverse[x, y, z, k]
                        num = np.float32(num + prior_increment)
                        if is_assigned_lesion:
                            num_lesion = np.float32(num_lesion + prior_increment)
                        energy = np.float32(energy + prior_objective[x, y, z, k])

                    energy_map[x, y, z, k] = energy
                    ordinary_reciprocal = (
                        np.float32(np.float32(1.0) / num)
                        if num > zero else infinity
                    )
                    if is_assigned_lesion:
                        lesion_reciprocal = (
                            np.float32(np.float32(1.0) / num_lesion)
                            if num_lesion > zero else infinity
                        )

                    direct = bool(np.any(closest[x, y, z, :2] == k))
                    neighbor = bool(np.any(closest[x, y, z, 2:] == k))
                    if kind == "LESION":
                        if is_assigned_lesion:
                            if direct:
                                factor = np.float32(1.0)
                            else:
                                factor = factors[x, y, z, k]
                            mixed = np.float32(
                                np.float32(np.float32(1.0) - lesion[x, y, z])
                                * ordinary_reciprocal
                            )
                            mixed = np.float32(
                                mixed + np.float32(lesion[x, y, z] * lesion_reciprocal)
                            )
                            current[x, y, z, k] = np.float32(factor * mixed)
                            lesion[x, y, z] = np.float32(
                                np.float32(factor * lesion[x, y, z]) * lesion_reciprocal
                            )
                        else:
                            # Selected Java quirk: the nonmatching NEIGHBOR
                            # route writes ordinary reciprocal without factor.
                            if direct or neighbor:
                                current[x, y, z, k] = ordinary_reciprocal
                            else:
                                current[x, y, z, k] = np.float32(
                                    factors[x, y, z, k] * ordinary_reciprocal
                                )
                            lesion[x, y, z] = np.float32(0.0)
                    else:
                        current[x, y, z, k] = np.float32(
                            factors[x, y, z, k] * ordinary_reciprocal
                        )
                    denominator = np.float32(denominator + current[x, y, z, k])
                    objective_at_voxel[k] = energy

                denominator_map[x, y, z] = denominator
                if denominator > np.float32(0.0):
                    for k, kind in enumerate(types):
                        current[x, y, z, k] = np.float32(current[x, y, z, k] / denominator)
                        if kind == "LESION" and labels[x, y, z] == k:
                            lesion[x, y, z] = np.float32(lesion[x, y, z] / denominator)
                        contribution = np.float32(
                            objective_at_voxel[k] * current[x, y, z, k]
                        )
                        contribution = np.float32(contribution * current[x, y, z, k])
                        objective = np.float32(objective + contribution)

    return ToadsLesionMembershipPass(
        class_data_terms=object_terms,
        lesion_data_terms=lesion_terms,
        class_energy=energy_map,
        memberships=current,
        lesion_membership=lesion,
        denominator=denominator_map,
        objective=objective,
    )
