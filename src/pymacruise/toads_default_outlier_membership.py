"""Fixed-input composition of the default TOADS atlas with one OUTLIER class."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES
from .toads_math import (
    toads_inverse_shape_prior_energy,
    toads_mask_intensity_energy,
    toads_optimized_tissue_data_term,
    toads_outlier_intensity_energy,
    toads_outlier_membership_mix,
    toads_shape_relation_factors,
)


_OPTIMIZED_TYPES = {"OPTIMIZED_CSF", "OPTIMIZED_GM", "OPTIMIZED_WM"}


@dataclass(frozen=True)
class ToadsDefaultOutlierMembershipPass:
    """Class/auxiliary data terms and one normalized membership pass."""

    class_data_energy: NDArray[np.float32]
    outlier_data_energy: NDArray[np.float32]
    memberships: NDArray[np.float32]
    outlier_membership: NDArray[np.float32]
    denominator: NDArray[np.float32]
    objective: np.float32


def toads_default_outlier_membership_pass(
    image_values: ArrayLike,
    initial_memberships: ArrayLike,
    initial_outlier_membership: ArrayLike,
    class_types: Sequence[str],
    centroids: ArrayLike,
    exponent_by_channel_class: ArrayLike,
    modality_weights: ArrayLike,
    class_variance_weights: ArrayLike,
    image_variances: ArrayLike,
    outlier_value: float,
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
) -> ToadsDefaultOutlierMembershipPass:
    """Compose selected default MASK/optimized classes and one OUTLIER.

    Removing the single OUTLIER name from ``class_types`` must leave the exact
    selected 11-class default atlas order. OUTLIER may occupy any slot; that
    slot's class centroid and weights are explicit inputs. The source's
    dangling-else behavior is retained: smoothing updates ordinary ``num``
    but not auxiliary ``numO``; the separate prior updates both. Statistics,
    initialization, topology and full iteration are out of scope. Spatial
    inputs must already be aligned in true RAS voxel order.
    """

    types = tuple(str(value).upper() for value in class_types)
    if types.count("OUTLIER") != 1 or tuple(
        kind for kind in types if kind != "OUTLIER"
    ) != DEFAULT_ATLAS_CLASS_TYPES:
        raise ValueError(
            "class_types must be the selected default 11-class order with one OUTLIER inserted"
        )

    images = np.asarray(image_values, dtype=np.float32)
    initial = np.asarray(initial_memberships, dtype=np.float32)
    initial_outlier = np.asarray(initial_outlier_membership, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    exponents = np.asarray(exponent_by_channel_class, dtype=np.float32)
    modality = np.asarray(modality_weights, dtype=np.float32)
    variance_weights = np.asarray(class_variance_weights, dtype=np.float32)
    image_variance = np.asarray(image_variances, dtype=np.float32)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    labels = np.asarray(closest_labels)
    smooth_prior = np.asarray(smoothing_priors, dtype=np.float32)
    similarity = np.asarray(centroid_similarity, dtype=np.float32)
    shape_prior = np.asarray(atlas_shape_priors, dtype=np.float32)

    if images.ndim != 4 or images.shape[-1] != 2 or min(images.shape[:3]) < 3:
        raise ValueError("image_values must be finite with shape (x>=3,y>=3,z>=3,2)")
    if not np.isfinite(images).all():
        raise ValueError("image_values must be finite")
    spatial, channels, classes = images.shape[:3], images.shape[-1], len(types)
    if initial.shape != spatial + (classes,) or not np.isfinite(initial).all():
        raise ValueError("initial_memberships must be finite and match the ordered class grid")
    if initial_outlier.shape != spatial or not np.isfinite(initial_outlier).all():
        raise ValueError("initial_outlier_membership must be finite and match the spatial grid")
    if centers.shape != (channels, classes) or not np.isfinite(centers).all():
        raise ValueError("centroids must be finite and match (channel,class) dimensions")
    if exponents.shape != centers.shape or not np.isfinite(exponents).all():
        raise ValueError("exponent_by_channel_class must be finite and match centroids")
    if variance_weights.shape != centers.shape or not np.isfinite(variance_weights).all():
        raise ValueError("class_variance_weights must be finite and match centroids")
    if modality.shape != (channels,) or image_variance.shape != (channels,):
        raise ValueError("modality_weights and image_variances must match channel count")
    if not np.isfinite(modality).all() or not np.isfinite(image_variance).all():
        raise ValueError("modality and image variance weights must be finite")
    if shapes.shape != initial.shape or not np.isfinite(shapes).all():
        raise ValueError("transformed_shapes must be finite and match memberships")
    if labels.shape != spatial + (4,) or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("closest_labels must be an integer four-slot map")
    if np.any((labels < -1) | (labels >= classes)):
        raise ValueError("closest_labels values must be -1 or valid ordered class indices")
    if smooth_prior.shape != (classes,) or not np.isfinite(smooth_prior).all():
        raise ValueError("smoothing_priors must be finite and match class count")
    if similarity.shape != (classes, classes) or not np.isfinite(similarity).all():
        raise ValueError("centroid_similarity must be finite and square by class")
    if shape_prior.shape != (classes,) or not np.isfinite(shape_prior).all():
        raise ValueError("atlas_shape_priors must be finite and match class count")
    coefficient = np.float32(prior_coefficient)
    smooth = np.float32(smoothing)
    factor = np.float32(prior_factor)
    outlier_scalar = np.float32(outlier_value)
    if not np.isfinite(coefficient) or not np.isfinite(smooth) or smooth < 0:
        raise ValueError("prior_coefficient and nonnegative smoothing must be finite")
    if not np.isfinite(factor) or not np.isfinite(outlier_scalar):
        raise ValueError("prior_factor and outlier_value must be finite")

    field = None if intensity_field is None else np.asarray(intensity_field, dtype=np.float32)
    bias = None if psf_bias is None else np.asarray(psf_bias, dtype=np.float32)
    if field is not None and (field.shape != images.shape or not np.isfinite(field).all()):
        raise ValueError("intensity_field must be finite and match image_values")
    if bias is not None and (bias.shape != images.shape or not np.isfinite(bias).all()):
        raise ValueError("psf_bias must be finite and match image_values")

    outlier_index = types.index("OUTLIER")
    mask_index = types.index("MASK")
    mask_energy = toads_mask_intensity_energy(
        images,
        modality,
        variance_weights[:, mask_index : mask_index + 1],
        image_variance,
    )
    ordinary_outlier, auxiliary_by_class = toads_outlier_intensity_energy(
        images,
        centers[:, outlier_index : outlier_index + 1],
        modality,
        variance_weights[:, outlier_index : outlier_index + 1],
        image_variance,
        float(outlier_scalar),
        intensity_field=field,
        psf_bias=bias,
    )
    class_energy = np.zeros(spatial + (classes,), dtype=np.float32)
    auxiliary_energy = np.zeros(spatial, dtype=np.float32)
    interior = (slice(1, -1), slice(1, -1), slice(1, -1))
    for index, kind in enumerate(types):
        if kind == "MASK":
            class_energy[interior + (index,)] = mask_energy[interior + (0,)]
        elif kind == "OUTLIER":
            auxiliary_energy[interior] = auxiliary_by_class[interior + (0,)]
            raw = ordinary_outlier[interior + (0,)]
            outlier_map_old = initial_outlier[interior]
            complement = np.asarray(np.float32(1.0) - outlier_map_old, dtype=np.float32)
            class_energy[interior + (index,)] = np.asarray(
                np.asarray(complement * raw, dtype=np.float32)
                + np.asarray(outlier_map_old * auxiliary_energy[interior], dtype=np.float32),
                dtype=np.float32,
            )
        else:
            selected = toads_optimized_tissue_data_term(
                images,
                centers[:, index : index + 1],
                exponents[:, index : index + 1],
                modality,
                variance_weights[:, index : index + 1],
                image_variance,
                intensity_field=field,
                psf_bias=bias,
            )
            class_energy[interior + (index,)] = selected[interior + (0,)]

    if not np.isfinite(class_energy).all() or not np.isfinite(auxiliary_energy).all():
        raise ValueError("class and auxiliary data-energy maps must be finite")

    if coefficient > np.float32(0.0):
        prior_inverse, prior_objective = toads_inverse_shape_prior_energy(
            shapes, shape_prior, similarity, coefficient
        )
        relation_shapes = shapes
    else:
        prior_inverse = np.zeros_like(class_energy)
        prior_objective = np.zeros_like(class_energy)
        relation_shapes = np.zeros_like(shapes)
    factors, suppressed = toads_shape_relation_factors(relation_shapes, labels, factor)
    class_energy[suppressed] = np.float32(0.0)
    auxiliary_energy[suppressed[..., outlier_index]] = np.float32(0.0)

    updated = initial.copy()
    outlier_map = initial_outlier.copy()
    denominator_map = np.zeros(spatial, dtype=np.float32)
    objective = np.float32(0.0)
    threshold, reciprocal_cap = np.float32(1e-10), np.float32(1e10)
    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    nx, ny, nz = spatial
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                denominator = np.float32(0.0)
                objective_at = np.zeros(classes, dtype=np.float32)
                for index, kind in enumerate(types):
                    if suppressed[x, y, z, index]:
                        updated[x, y, z, index] = np.float32(0.0)
                        continue
                    num = class_energy[x, y, z, index]
                    num_outlier = auxiliary_energy[x, y, z]
                    current_outlier = outlier_map[x, y, z]
                    objective_value = num
                    if smooth > np.float32(0.0):
                        neighbor_sum = np.float32(0.0)
                        for dx, dy, dz in offsets:
                            for other in range(classes):
                                if other == index:
                                    continue
                                neighbor = updated[x + dx, y + dy, z + dz, other]
                                neighbor_sum = np.float32(
                                    neighbor_sum + np.float32(neighbor * neighbor)
                                )
                        smoothing_increment = np.float32(smooth_prior[index] * smooth)
                        smoothing_increment = np.float32(smoothing_increment * neighbor_sum)
                        num = np.float32(num + smoothing_increment)
                        objective_increment = np.float32(
                            np.float32(np.float32(0.5) * smooth_prior[index]) * smooth
                        )
                        objective_increment = np.float32(objective_increment * neighbor_sum)
                        objective_value = np.float32(objective_value + objective_increment)
                    num = np.float32(num + prior_inverse[x, y, z, index])
                    if kind == "OUTLIER":
                        num_outlier = np.float32(
                            num_outlier + prior_inverse[x, y, z, index]
                        )
                    objective_value = np.float32(
                        objective_value + prior_objective[x, y, z, index]
                    )
                    reciprocal = (
                        np.float32(np.float32(1.0) / num)
                        if num > threshold else reciprocal_cap
                    )
                    if kind == "OUTLIER":
                        reciprocal_outlier = (
                            np.float32(np.float32(1.0) / num_outlier)
                            if num_outlier > threshold else reciprocal_cap
                        )
                        membership, outlier_map[x, y, z] = toads_outlier_membership_mix(
                            reciprocal,
                            current_outlier,
                            reciprocal_outlier,
                            factors[x, y, z, index],
                        )
                    else:
                        membership = np.float32(factors[x, y, z, index] * reciprocal)
                    updated[x, y, z, index] = membership
                    denominator = np.float32(denominator + membership)
                    objective_at[index] = objective_value

                denominator_map[x, y, z] = denominator
                if denominator > np.float32(0.0):
                    for index in range(classes):
                        updated[x, y, z, index] = np.float32(
                            updated[x, y, z, index] / denominator
                        )
                        if index == outlier_index:
                            outlier_map[x, y, z] = np.float32(
                                outlier_map[x, y, z] / denominator
                            )
                        contribution = np.float32(objective_at[index] * updated[x, y, z, index])
                        contribution = np.float32(contribution * updated[x, y, z, index])
                        objective = np.float32(objective + contribution)

    return ToadsDefaultOutlierMembershipPass(
        class_energy,
        auxiliary_energy,
        updated,
        outlier_map,
        denominator_map,
        objective,
    )
