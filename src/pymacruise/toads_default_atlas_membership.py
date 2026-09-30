"""Fixed-input exact-membership pass for the selected 11-class TOADS atlas.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``. This
composes one fixed-statistics tissue-classification stage; it does not perform
the papers' complete alternating topology/registration workflow.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_math import (
    toads_inverse_shape_prior_energy,
    toads_mask_intensity_energy,
    toads_optimized_tissue_data_term,
    toads_shape_relation_factors,
)


DEFAULT_ATLAS_CLASS_TYPES = (
    "MASK",
    "OPTIMIZED_CSF", "OPTIMIZED_CSF",
    "OPTIMIZED_GM", "OPTIMIZED_GM", "OPTIMIZED_GM", "OPTIMIZED_GM", "OPTIMIZED_GM",
    "OPTIMIZED_WM", "OPTIMIZED_WM", "OPTIMIZED_WM",
)
_OPTIMIZED_TYPES = {"OPTIMIZED_CSF", "OPTIMIZED_GM", "OPTIMIZED_WM"}


@dataclass(frozen=True)
class ToadsDefaultAtlasMembershipPass:
    """Class data terms, updated maps, denominator, and objective."""

    class_data_energy: NDArray[np.float32]
    memberships: NDArray[np.float32]
    denominator: NDArray[np.float32]
    objective: np.float32


def toads_default_atlas_membership_pass(
    image_values: ArrayLike,
    initial_memberships: ArrayLike,
    centroids: ArrayLike,
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
    class_types: Sequence[str] = DEFAULT_ATLAS_CLASS_TYPES,
    intensity_field: ArrayLike | None = None,
    psf_bias: ArrayLike | None = None,
) -> ToadsDefaultAtlasMembershipPass:
    """Compose the selected default-atlas pass with fixed statistics.

    The only accepted class order is the 11-class descriptor order used by
    ``MedicAlgorithmToads08``: one MASK, two optimized CSF, five optimized GM,
    and three optimized WM classes. Repeated optimized classes have distinct
    centroids and variance weights but share the source channel/category
    PowerTable exponents. Statistics, transformed shapes, relation labels,
    similarity, and atlas priors are explicit inputs; no TOADS estimation,
    topology or outer iteration is included.

    Source `computeExactMemberships()` leaves its local shape vector zero when
    the prior coefficient is not positive. Thus prior-only classes are
    suppressed in that mode even when transformed-shape inputs were supplied.
    This adapter preserves that behavior and the in-place x/y/z/class sweep.
    Inputs are already in true-RAS voxel order.
    """

    types = tuple(str(value).upper() for value in class_types)
    if types != DEFAULT_ATLAS_CLASS_TYPES:
        raise ValueError("class_types must match the selected default 11-class atlas order")
    images = np.asarray(image_values, dtype=np.float32)
    initial = np.asarray(initial_memberships, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    exponents = np.asarray(exponent_by_channel_class, dtype=np.float32)
    variances = np.asarray(class_variance_weights, dtype=np.float32)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    labels = np.asarray(closest_labels)
    if images.ndim != 4 or images.shape[-1] != 2 or min(images.shape[:3]) < 3:
        raise ValueError("image_values must have shape (x>=3,y>=3,z>=3,2)")
    spatial = images.shape[:3]
    classes = len(types)
    if initial.shape != spatial + (classes,) or not np.isfinite(initial).all():
        raise ValueError("initial_memberships must match the 11-class image grid")
    if centers.shape != (2, classes) or not np.isfinite(centers).all():
        raise ValueError("centroids must be finite with shape (2,11)")
    if exponents.shape != centers.shape or not np.isfinite(exponents).all():
        raise ValueError("exponent_by_channel_class must be finite with shape (2,11)")
    if variances.shape != centers.shape or not np.isfinite(variances).all():
        raise ValueError("class_variance_weights must be finite with shape (2,11)")
    if shapes.shape != initial.shape or not np.isfinite(shapes).all():
        raise ValueError("transformed_shapes must be finite and match initial_memberships")
    if labels.shape != spatial + (4,) or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("closest_labels must be an integer four-slot map")
    modality = np.asarray(modality_weights, dtype=np.float32)
    image_var = np.asarray(image_variances, dtype=np.float32)
    smooth_prior = np.asarray(smoothing_priors, dtype=np.float32)
    similarity = np.asarray(centroid_similarity, dtype=np.float32)
    shape_prior = np.asarray(atlas_shape_priors, dtype=np.float32)
    if modality.shape != (2,) or image_var.shape != (2,):
        raise ValueError("modality_weights and image_variances must have shape (2,)")
    if smooth_prior.shape != (classes,) or not np.isfinite(smooth_prior).all():
        raise ValueError("smoothing_priors must have 11 finite values")
    if similarity.shape != (classes, classes) or not np.isfinite(similarity).all():
        raise ValueError("centroid_similarity must be finite and 11 by 11")
    if shape_prior.shape != (classes,) or not np.isfinite(shape_prior).all():
        raise ValueError("atlas_shape_priors must contain 11 finite values")
    if not np.isfinite(images).all() or not np.isfinite(modality).all() or not np.isfinite(image_var).all():
        raise ValueError("image and intensity weights must be finite")
    for kind in _OPTIMIZED_TYPES:
        columns = [i for i, value in enumerate(types) if value == kind]
        if len(columns) > 1 and not np.array_equal(
            exponents[:, columns], np.broadcast_to(exponents[:, columns[0], None], (2, len(columns)))
        ):
            raise ValueError(f"all {kind} classes must share category/channel exponents")

    field = None if intensity_field is None else np.asarray(intensity_field, dtype=np.float32)
    bias = None if psf_bias is None else np.asarray(psf_bias, dtype=np.float32)
    if field is not None and (field.shape != images.shape or not np.isfinite(field).all()):
        raise ValueError("intensity_field must be finite and match image_values")
    if bias is not None and (bias.shape != images.shape or not np.isfinite(bias).all()):
        raise ValueError("psf_bias must be finite and match image_values")

    class_energy = np.zeros(spatial + (classes,), dtype=np.float32)
    mask_energy = toads_mask_intensity_energy(images, modality, variances, image_var)
    class_energy[1:-1, 1:-1, 1:-1, 0] = mask_energy[1:-1, 1:-1, 1:-1, 0]
    for class_index, kind in enumerate(types):
        if kind not in _OPTIMIZED_TYPES:
            continue
        category_energy = toads_optimized_tissue_data_term(
            images,
            centers[:, class_index : class_index + 1],
            exponents[:, class_index : class_index + 1],
            modality,
            variances[:, class_index : class_index + 1],
            image_var,
            intensity_field=field,
            psf_bias=bias,
        )
        class_energy[1:-1, 1:-1, 1:-1, class_index] = category_energy[
            1:-1, 1:-1, 1:-1, 0
        ]

    coefficient = np.float32(prior_coefficient)
    if not np.isfinite(coefficient):
        raise ValueError("prior_coefficient must be finite")
    if coefficient > np.float32(0.0):
        inverse_prior, objective_prior = toads_inverse_shape_prior_energy(
            shapes, shape_prior, similarity, coefficient
        )
        relation_shapes = shapes
    else:
        inverse_prior = np.zeros_like(class_energy)
        objective_prior = np.zeros_like(class_energy)
        relation_shapes = np.zeros_like(shapes)

    factors, suppressed = toads_shape_relation_factors(
        relation_shapes, labels, prior_factor
    )
    class_energy[suppressed] = np.float32(0.0)

    smooth = np.float32(smoothing)
    if not np.isfinite(smooth) or smooth < np.float32(0.0):
        raise ValueError("smoothing must be finite and nonnegative")
    prior_scale = np.float32(prior_factor)
    if not np.isfinite(prior_scale):
        raise ValueError("prior_factor must be finite")

    updated = initial.copy()
    denominator_map = np.zeros(spatial, dtype=np.float32)
    objective = np.float32(0.0)
    threshold, reciprocal_cap = np.float32(1.0e-10), np.float32(1.0e10)
    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    nx, ny, nz = spatial
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                denominator = np.float32(0.0)
                objective_at = np.zeros(classes, dtype=np.float32)
                for class_index in range(classes):
                    # computeExactMemberships() skips data and updates for
                    # prior-only classes with zero shape.
                    if suppressed[x, y, z, class_index]:
                        updated[x, y, z, class_index] = np.float32(0.0)
                        continue
                    num = class_energy[x, y, z, class_index]
                    objective_value = num
                    if smooth > np.float32(0.0):
                        neighbor_sum = np.float32(0.0)
                        for dx, dy, dz in offsets:
                            for other in range(classes):
                                if other == class_index:
                                    continue
                                value = updated[x + dx, y + dy, z + dz, other]
                                neighbor_sum = np.float32(
                                    neighbor_sum + np.float32(value * value)
                                )
                        increment = np.float32(
                            np.float32(smooth_prior[class_index] * smooth) * neighbor_sum
                        )
                        num = np.float32(num + increment)
                        objective_increment = np.float32(
                            np.float32(np.float32(0.5) * smooth_prior[class_index]) * smooth
                        )
                        objective_increment = np.float32(objective_increment * neighbor_sum)
                        objective_value = np.float32(objective_value + objective_increment)
                    num = np.float32(num + inverse_prior[x, y, z, class_index])
                    objective_value = np.float32(
                        objective_value + objective_prior[x, y, z, class_index]
                    )
                    reciprocal = (
                        np.float32(np.float32(1.0) / num)
                        if num > threshold else reciprocal_cap
                    )
                    value = np.float32(factors[x, y, z, class_index] * reciprocal)
                    updated[x, y, z, class_index] = value
                    denominator = np.float32(denominator + value)
                    objective_at[class_index] = objective_value
                denominator_map[x, y, z] = denominator
                if denominator > np.float32(0.0):
                    for class_index in range(classes):
                        updated[x, y, z, class_index] = np.float32(
                            updated[x, y, z, class_index] / denominator
                        )
                        contribution = np.float32(
                            objective_at[class_index] * updated[x, y, z, class_index]
                        )
                        contribution = np.float32(
                            contribution * updated[x, y, z, class_index]
                        )
                        objective = np.float32(objective + contribution)

    return ToadsDefaultAtlasMembershipPass(class_energy, updated, denominator_map, objective)
