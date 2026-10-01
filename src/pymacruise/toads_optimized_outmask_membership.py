"""Fixed-input base-TOADS optimized-tissue plus OUTMASK membership pass.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``, describe
statistical tissue classification with atlas-informed constraints. This
fixed-input pass composes selected optimized-tissue classes and OUTMASK in
that broad role. Its source-specific energy terms, order, and masking
semantics come from VandyCRUISE, not directly from the papers; this is not a
complete TOADS iteration.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_math import (
    toads_inverse_shape_prior_energy,
    toads_mask_intensity_energy,
    toads_object_intensity_energy,
    toads_optimized_tissue_data_term,
    toads_outlier_intensity_energy,
    toads_outlier_membership_mix,
    toads_shape_relation_factors,
)


_STANDARD_TYPES = {"OBJECT", "MASK"}
_OPTIMIZED_TYPES = {"OPTIMIZED_CSF", "OPTIMIZED_GM", "OPTIMIZED_WM"}
_SUPPORTED_TYPES = _STANDARD_TYPES | _OPTIMIZED_TYPES | {"OUTMASK"}


@dataclass(frozen=True)
class ToadsOptimizedOutmaskMembershipPass:
    """Energy maps, normalized class/auxiliary maps, and objective."""

    class_data_energy: NDArray[np.float32]
    outlier_data_energy: NDArray[np.float32]
    memberships: NDArray[np.float32]
    outlier_membership: NDArray[np.float32]
    denominator: NDArray[np.float32]
    objective: np.float32


def toads_optimized_outmask_membership_pass(
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
) -> ToadsOptimizedOutmaskMembershipPass:
    """Compose OBJECT/MASK/optimized tissues with exactly one OUTMASK.

    Statistics and maps are explicit fixed inputs. Optimized classes use their
    selected category exponent in the already source-matched PowerTable data
    kernel. The OUTMASK class's ordinary term is mask-distance energy and its
    auxiliary term is the fixed outlier energy mixture. The selected Java
    unbraced smoothing conditional adds smoothing to ordinary class energy
    but skips auxiliary OUTMASK smoothing; its separate prior branch updates
    both. Membership and objective updates preserve class order and the
    selected six-neighbor in-place sweep.

    Spatial arrays are already aligned in RAS voxel order. This composes only
    this fixed-input family; LESION, OUTLIER, RAYLEIGH, estimation, topology,
    and outer iteration are not implemented here.
    """

    types = tuple(str(value).upper() for value in class_types)
    if (
        not types
        or types.count("OUTMASK") != 1
        or not set(types).intersection(_STANDARD_TYPES)
        or not set(types).intersection(_OPTIMIZED_TYPES)
        or any(kind not in _SUPPORTED_TYPES for kind in types)
    ):
        raise ValueError(
            "class_types must contain OBJECT or MASK, at least one optimized tissue, "
            "and exactly one OUTMASK"
        )

    images = np.asarray(image_values, dtype=np.float32)
    initial = np.asarray(initial_memberships, dtype=np.float32)
    initial_outlier = np.asarray(initial_outlier_membership, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    exponents = np.asarray(exponent_by_channel_class, dtype=np.float32)
    variance_weights = np.asarray(class_variance_weights, dtype=np.float32)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    labels = np.asarray(closest_labels)
    if images.ndim != 4 or images.shape[-1] < 1 or not np.isfinite(images).all():
        raise ValueError("image_values must be finite with shape (x,y,z,modality)")
    spatial = images.shape[:3]
    channels, classes = images.shape[-1], len(types)
    if initial.shape != spatial + (classes,):
        raise ValueError("initial_memberships must match the class grid")
    if initial_outlier.shape != spatial:
        raise ValueError("initial_outlier_membership must match the spatial grid")
    if centers.shape != (channels, classes):
        raise ValueError("centroids must match (modality,class) dimensions")
    if exponents.shape != centers.shape:
        raise ValueError("exponent_by_channel_class must match centroid dimensions")
    if variance_weights.shape != centers.shape:
        raise ValueError("class_variance_weights must match centroid dimensions")
    if shapes.shape != initial.shape:
        raise ValueError("transformed_shapes must match initial_memberships")
    if labels.shape != spatial + (4,):
        raise ValueError("closest_labels must match the spatial grid with four slots")
    if not np.isfinite(initial).all() or not np.isfinite(initial_outlier).all():
        raise ValueError("initial membership maps must be finite")
    if not np.isfinite(shapes).all() or np.any(shapes[1:-1, 1:-1, 1:-1] == 0):
        raise ValueError("all class shapes must be finite and active in the interior")
    outlier_value32 = np.float32(outlier_value)
    if not np.isfinite(outlier_value32):
        raise ValueError("outlier_value must be finite")

    object_energy = toads_object_intensity_energy(
        images, centers, modality_weights, variance_weights, image_variances,
        intensity_field=intensity_field, psf_bias=psf_bias,
    )
    mask_energy = toads_mask_intensity_energy(
        images, modality_weights, variance_weights, image_variances,
    )
    _, outlier_by_class = toads_outlier_intensity_energy(
        images, centers, modality_weights, variance_weights, image_variances,
        float(outlier_value32), intensity_field=intensity_field, psf_bias=psf_bias,
    )
    class_energy = np.zeros(spatial + (classes,), dtype=np.float32)
    auxiliary_energy = np.zeros(spatial, dtype=np.float32)
    interior = (slice(1, -1), slice(1, -1), slice(1, -1))
    outmask_index = types.index("OUTMASK")
    for index, kind in enumerate(types):
        if kind == "OBJECT":
            class_energy[interior + (index,)] = object_energy[interior + (index,)]
        elif kind == "MASK" or kind == "OUTMASK":
            class_energy[interior + (index,)] = mask_energy[interior + (index,)]
        else:
            class_power = toads_optimized_tissue_data_term(
                images,
                centers[:, index : index + 1],
                exponents[:, index : index + 1],
                modality_weights,
                variance_weights[:, index : index + 1],
                image_variances,
                intensity_field=intensity_field,
                psf_bias=psf_bias,
            )
            class_energy[interior + (index,)] = class_power[interior + (0,)]
    auxiliary_energy[interior] = outlier_by_class[interior + (outmask_index,)]

    prior_inverse, prior_objective = toads_inverse_shape_prior_energy(
        shapes, atlas_shape_priors, centroid_similarity, prior_coefficient
    )
    factors, suppressed = toads_shape_relation_factors(shapes, labels, prior_factor)
    if np.any(suppressed[1:-1, 1:-1, 1:-1]):
        raise ValueError("all class shapes must be active in the interior")
    if not np.isfinite(class_energy).all() or not np.isfinite(auxiliary_energy).all():
        raise ValueError("class and auxiliary energies must be finite")

    updated = initial.copy()
    outlier = initial_outlier.copy()
    smooth = np.float32(smoothing)
    zero = np.float32(1e-10)
    infinity = np.float32(1e10)
    if not np.isfinite(smooth) or smooth < 0:
        raise ValueError("smoothing must be finite and nonnegative")
    total_objective = np.float32(0.0)
    denominators = np.zeros(spatial, dtype=np.float32)
    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    nx, ny, nz = spatial
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                denominator = np.float32(0.0)
                objective_at = np.zeros(classes, dtype=np.float32)
                for index, kind in enumerate(types):
                    num = class_energy[x, y, z, index]
                    num_outlier = auxiliary_energy[x, y, z]
                    objective_value = num
                    old_outlier = outlier[x, y, z]
                    if kind == "OUTMASK":
                        objective_value = np.float32(
                            np.float32(np.float32(1.0) - old_outlier) * num
                            + old_outlier * num_outlier
                        )

                    if smooth > np.float32(0.0):
                        neighbor_sum = np.float32(0.0)
                        for dx, dy, dz in offsets:
                            for other in range(classes):
                                if other == index:
                                    continue
                                value = updated[x + dx, y + dy, z + dz, other]
                                neighbor_sum = np.float32(
                                    neighbor_sum + np.float32(value * value)
                                )
                        increment = np.float32(np.float32(smoothing_priors[index]) * smooth)
                        increment = np.float32(increment * neighbor_sum)
                        num = np.float32(num + increment)
                        objective_value = np.float32(
                            objective_value
                            + np.float32(
                                np.float32(
                                    np.float32(0.5) * np.float32(smoothing_priors[index])
                                )
                                * smooth
                            )
                            * neighbor_sum
                        )
                        # The selected dangling-else skips OUTMASK auxiliary
                        # smoothing, even though `num` above receives smoothing.

                    if np.float32(prior_coefficient) > np.float32(0.0):
                        num = np.float32(num + prior_inverse[x, y, z, index])
                        if kind == "OUTMASK":
                            num_outlier = np.float32(
                                num_outlier + prior_inverse[x, y, z, index]
                            )
                        objective_value = np.float32(
                            objective_value + prior_objective[x, y, z, index]
                        )

                    reciprocal = np.float32(1.0) / num if num > zero else infinity
                    if kind == "OUTMASK":
                        reciprocal_outlier = (
                            np.float32(1.0) / num_outlier if num_outlier > zero else infinity
                        )
                        value, outlier[x, y, z] = toads_outlier_membership_mix(
                            reciprocal, old_outlier, reciprocal_outlier, factors[x, y, z, index]
                        )
                    else:
                        value = np.float32(factors[x, y, z, index] * reciprocal)
                    updated[x, y, z, index] = value
                    denominator = np.float32(denominator + value)
                    objective_at[index] = objective_value

                if denominator > np.float32(0.0):
                    denominators[x, y, z] = denominator
                    for index in range(classes):
                        updated[x, y, z, index] = np.float32(
                            updated[x, y, z, index] / denominator
                        )
                        if index == outmask_index:
                            outlier[x, y, z] = np.float32(outlier[x, y, z] / denominator)
                        contribution = np.float32(
                            objective_at[index] * updated[x, y, z, index]
                        )
                        contribution = np.float32(contribution * updated[x, y, z, index])
                        total_objective = np.float32(total_objective + contribution)

    return ToadsOptimizedOutmaskMembershipPass(
        class_energy, auxiliary_energy, updated, outlier, denominators, total_objective
    )
