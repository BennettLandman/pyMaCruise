"""Fixed-input membership pass for the selected lesion-bearing TOADS atlas."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_atlas_type_membership import (
    map_toads_topologies,
    read_toads_atlas_class_types,
)
from .toads_math import (
    toads_inverse_shape_prior_energy,
    toads_mask_intensity_energy,
    toads_object_intensity_energy,
    toads_optimized_tissue_data_term,
    toads_shape_relation_factors,
)


LESION_TOADS_CLASS_TYPES = (
    "MASK", "OPTIMIZED_CSF", "OPTIMIZED_CSF",
    "OBJECT", "OBJECT", "OBJECT", "OBJECT", "OBJECT", "OBJECT", "OBJECT",
    "LESION", "OBJECT", "OBJECT",
)
LESION_CLASS_INDEX = 10
_OPTIMIZED_TYPES = {"OPTIMIZED_CSF", "OPTIMIZED_GM", "OPTIMIZED_WM"}


@dataclass(frozen=True)
class ToadsLesionAtlasTypeMembershipPass:
    """Selected class/data terms and one updated exact-membership pass."""

    class_types: tuple[str, ...]
    class_data_terms: NDArray[np.float32]
    lesion_data_terms: NDArray[np.float32]
    memberships: NDArray[np.float32]
    lesion_membership: NDArray[np.float32]
    denominator: NDArray[np.float32]
    objective: np.float32


def read_lesion_toads_atlas_class_types(path: str | Path) -> tuple[str, ...]:
    """Read descriptor order and require the selected lesion 13-class mix."""
    types = read_toads_atlas_class_types(path)
    if types != LESION_TOADS_CLASS_TYPES:
        raise ValueError("descriptor must match selected 13-class LesionTOADS topology order")
    return types


def map_lesion_toads_topologies(topologies: Sequence[str]) -> tuple[str, ...]:
    """Apply selected case-sensitive topology mapping and check its order."""
    types = map_toads_topologies(topologies)
    if types != LESION_TOADS_CLASS_TYPES:
        raise ValueError("topologies must match selected 13-class LesionTOADS order")
    return types


def toads_lesion_atlas_type_membership_pass(
    image_values: ArrayLike,
    initial_memberships: ArrayLike,
    initial_lesion_membership: ArrayLike,
    segmentation: ArrayLike,
    class_types: Sequence[str],
    centroids: ArrayLike,
    lesion_centroid: ArrayLike,
    exponent_by_channel_class: ArrayLike,
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
) -> ToadsLesionAtlasTypeMembershipPass:
    """Run one source-ordered 13-class exact-membership pass with fixed state.

    The LESION category uses its ordinary centroid residual for every active
    voxel, and its second lesion-centroid term only when the hard segmentation
    equals class index 10. The lesion-membership map is updated in the same
    class slot and divided by the common membership denominator. Input arrays
    are true-RAS subject-grid arrays; atlas template sampling/statistics and
    all workflow control are outside this adapter.
    """
    types = tuple(str(value).upper() for value in class_types)
    if types != LESION_TOADS_CLASS_TYPES:
        raise ValueError("class_types must match selected 13-class LesionTOADS order")
    images = np.asarray(image_values, dtype=np.float32)
    initial = np.asarray(initial_memberships, dtype=np.float32)
    lesion_initial = np.asarray(initial_lesion_membership, dtype=np.float32)
    labels = np.asarray(segmentation)
    centers = np.asarray(centroids, dtype=np.float32)
    lesion_center = np.asarray(lesion_centroid, dtype=np.float32)
    exponents = np.asarray(exponent_by_channel_class, dtype=np.float32)
    variance_weights = np.asarray(class_variance_weights, dtype=np.float32)
    variances = np.asarray(image_variances, dtype=np.float32)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    closest = np.asarray(closest_labels)
    classes = len(types)
    if images.ndim != 4 or images.shape[-1] < 1 or min(images.shape[:3]) < 3:
        raise ValueError("image_values must have shape (x>=3,y>=3,z>=3,modalities>=1)")
    spatial, channels = images.shape[:3], images.shape[-1]
    if initial.shape != spatial + (classes,) or not np.isfinite(initial).all():
        raise ValueError("initial_memberships must be finite and match the 13-class grid")
    if lesion_initial.shape != spatial or not np.isfinite(lesion_initial).all():
        raise ValueError("initial_lesion_membership must match the subject grid")
    if labels.shape != spatial or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer map matching the subject grid")
    if centers.shape != (channels, classes) or not np.isfinite(centers).all():
        raise ValueError("centroids must be finite with shape (channels,13)")
    if lesion_center.shape != (channels,) or not np.isfinite(lesion_center).all():
        raise ValueError("lesion_centroid must contain one finite value per channel")
    if exponents.shape != (channels, classes) or not np.isfinite(exponents).all():
        raise ValueError("exponent_by_channel_class must have shape (channels,13)")
    if variance_weights.shape != (channels, classes) or not np.isfinite(variance_weights).all():
        raise ValueError("class_variance_weights must have shape (channels,13)")
    if variances.shape != (channels,) or not np.isfinite(variances).all():
        raise ValueError("image_variances must contain one finite value per channel")
    if shapes.shape != initial.shape or not np.isfinite(shapes).all():
        raise ValueError("transformed_shapes must be finite and match memberships")
    if closest.shape != spatial + (4,) or not np.issubdtype(closest.dtype, np.integer):
        raise ValueError("closest_labels must be an integer four-slot subject-grid map")
    if np.any((labels < -1) | (labels >= classes)):
        raise ValueError("segmentation labels must be -1 or valid class indices")
    weights = np.asarray(modality_weights, dtype=np.float32)
    if weights.shape != (channels, 2) or not np.isfinite(weights).all():
        raise ValueError("modality_weights must contain ordinary and lesion rows per channel")
    smooth_priors = np.asarray(smoothing_priors, dtype=np.float32)
    similarities = np.asarray(centroid_similarity, dtype=np.float32)
    shape_priors = np.asarray(atlas_shape_priors, dtype=np.float32)
    if smooth_priors.shape != (classes,) or not np.isfinite(smooth_priors).all():
        raise ValueError("smoothing_priors must contain 13 finite values")
    if similarities.shape != (classes, classes) or not np.isfinite(similarities).all():
        raise ValueError("centroid_similarity must be finite and 13 by 13")
    if shape_priors.shape != (classes,) or not np.isfinite(shape_priors).all():
        raise ValueError("atlas_shape_priors must contain 13 finite values")
    if not np.isfinite(images).all():
        raise ValueError("image_values must be finite")

    for kind in _OPTIMIZED_TYPES:
        columns = [i for i, value in enumerate(types) if value == kind]
        if len(columns) > 1 and not np.array_equal(
            exponents[:, columns],
            np.broadcast_to(exponents[:, columns[0], None], (channels, len(columns))),
        ):
            raise ValueError(f"all {kind} classes must share category/channel exponents")
    field = None if intensity_field is None else np.asarray(intensity_field, dtype=np.float32)
    bias = None if psf_bias is None else np.asarray(psf_bias, dtype=np.float32)
    if field is not None and (field.shape != images.shape or not np.isfinite(field).all()):
        raise ValueError("intensity_field must be finite and match image_values")
    if bias is not None and (bias.shape != images.shape or not np.isfinite(bias).all()):
        raise ValueError("psf_bias must be finite and match image_values")
    smooth = np.float32(smoothing)
    prior_scale = np.float32(prior_factor)
    coefficient = np.float32(prior_coefficient)
    if not np.isfinite([smooth, prior_scale, coefficient]).all() or smooth < 0:
        raise ValueError("membership coefficients must be finite and smoothing nonnegative")

    class_terms = np.zeros(spatial + (classes,), dtype=np.float32)
    lesion_terms = np.zeros_like(class_terms)
    interior = (slice(1, -1), slice(1, -1), slice(1, -1))
    corrected = images if field is None else np.asarray(field * images, dtype=np.float32)
    if bias is not None:
        corrected = np.asarray(corrected - bias, dtype=np.float32)
    for k, kind in enumerate(types):
        weights_ordinary = weights[:, 0]
        one_variance = variance_weights[:, k : k + 1]
        if kind == "MASK":
            energy = toads_mask_intensity_energy(images, weights_ordinary, one_variance, variances)
        elif kind == "OPTIMIZED_CSF":
            energy = toads_optimized_tissue_data_term(
                images, centers[:, k : k + 1], exponents[:, k : k + 1],
                weights_ordinary, one_variance, variances,
                intensity_field=field, psf_bias=bias,
            )
        else:
            energy = toads_object_intensity_energy(
                images, centers[:, k : k + 1], weights_ordinary,
                one_variance, variances, intensity_field=field, psf_bias=bias,
            )
        class_terms[interior + (k,)] = energy[interior + (0,)]
    lesion_index = LESION_CLASS_INDEX
    for x in range(1, spatial[0] - 1):
        for y in range(1, spatial[1] - 1):
            for z in range(1, spatial[2] - 1):
                if labels[x, y, z] != lesion_index:
                    continue
                value = np.float32(0.0)
                for c in range(channels):
                    residual = np.float32(corrected[x, y, z, c] - lesion_center[c])
                    weight = np.float32(weights[c, 1] * variance_weights[c, lesion_index])
                    weight = np.float32(weight * variances[c])
                    contribution = np.float32(weight * residual)
                    contribution = np.float32(contribution * residual)
                    value = np.float32(value + contribution)
                lesion_terms[x, y, z, lesion_index] = value

    relation_shapes = shapes if coefficient > 0 else np.zeros_like(shapes)
    factors, suppressed = toads_shape_relation_factors(relation_shapes, closest, prior_scale)
    class_terms[suppressed] = np.float32(0.0)
    lesion_terms[suppressed] = np.float32(0.0)
    inverse_prior, objective_prior = toads_inverse_shape_prior_energy(
        shapes, shape_priors, similarities, coefficient
    )
    current = initial.copy()
    lesion = lesion_initial.copy()
    denominator_map = np.zeros(spatial, dtype=np.float32)
    objective = np.float32(0.0)
    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    zero, infinity = np.float32(1e-10), np.float32(1e10)

    for x in range(1, spatial[0] - 1):
        for y in range(1, spatial[1] - 1):
            for z in range(1, spatial[2] - 1):
                denominator = np.float32(0.0)
                objective_at = np.zeros(classes, dtype=np.float32)
                for k, kind in enumerate(types):
                    if suppressed[x, y, z, k]:
                        current[x, y, z, k] = np.float32(0.0)
                        continue
                    assigned_lesion = kind == "LESION" and labels[x, y, z] == k
                    num = class_terms[x, y, z, k]
                    num_lesion = lesion_terms[x, y, z, k] if assigned_lesion else np.float32(0.0)
                    if assigned_lesion:
                        energy = np.float32(
                            np.float32(np.float32(1.0) - lesion[x, y, z]) * num
                            + np.float32(lesion[x, y, z] * num_lesion)
                        )
                    else:
                        energy = num
                    if smooth > 0:
                        neighbor_sum = np.float32(0.0)
                        for dx, dy, dz in offsets:
                            for other in range(classes):
                                if other != k:
                                    q = current[x + dx, y + dy, z + dz, other]
                                    neighbor_sum = np.float32(neighbor_sum + np.float32(q * q))
                        increment = np.float32(np.float32(smooth_priors[k] * smooth) * neighbor_sum)
                        num = np.float32(num + increment)
                        if assigned_lesion:
                            num_lesion = np.float32(num_lesion + increment)
                        energy = np.float32(
                            energy + np.float32(
                                np.float32(np.float32(np.float32(0.5) * smooth_priors[k]) * smooth)
                                * neighbor_sum
                            )
                        )
                    if coefficient > 0:
                        prior = inverse_prior[x, y, z, k]
                        num = np.float32(num + prior)
                        if assigned_lesion:
                            num_lesion = np.float32(num_lesion + prior)
                        energy = np.float32(energy + objective_prior[x, y, z, k])
                    ordinary_reciprocal = np.float32(1.0) / num if num > zero else infinity
                    if assigned_lesion:
                        lesion_reciprocal = np.float32(1.0) / num_lesion if num_lesion > zero else infinity
                    direct = bool(np.any(closest[x, y, z, :2] == k))
                    neighbor = bool(np.any(closest[x, y, z, 2:] == k))
                    if kind == "LESION":
                        if assigned_lesion:
                            factor = np.float32(1.0) if direct else factors[x, y, z, k]
                            value = np.float32(np.float32(np.float32(1.0) - lesion[x, y, z]) * ordinary_reciprocal)
                            value = np.float32(value + np.float32(lesion[x, y, z] * lesion_reciprocal))
                            current[x, y, z, k] = np.float32(factor * value)
                            lesion[x, y, z] = np.float32(np.float32(factor * lesion[x, y, z]) * lesion_reciprocal)
                        else:
                            current[x, y, z, k] = (
                                ordinary_reciprocal if direct or neighbor
                                else np.float32(factors[x, y, z, k] * ordinary_reciprocal)
                            )
                            lesion[x, y, z] = np.float32(0.0)
                    else:
                        current[x, y, z, k] = np.float32(factors[x, y, z, k] * ordinary_reciprocal)
                    denominator = np.float32(denominator + current[x, y, z, k])
                    objective_at[k] = energy
                denominator_map[x, y, z] = denominator
                if denominator > 0:
                    for k, kind in enumerate(types):
                        current[x, y, z, k] = np.float32(current[x, y, z, k] / denominator)
                        if kind == "LESION" and labels[x, y, z] == k:
                            lesion[x, y, z] = np.float32(lesion[x, y, z] / denominator)
                        contribution = np.float32(objective_at[k] * current[x, y, z, k])
                        contribution = np.float32(contribution * current[x, y, z, k])
                        objective = np.float32(objective + contribution)

    return ToadsLesionAtlasTypeMembershipPass(
        types, class_terms, lesion_terms, current, lesion, denominator_map, objective
    )
