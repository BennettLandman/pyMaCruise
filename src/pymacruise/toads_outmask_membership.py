"""Fixed-input base-TOADS OBJECT/MASK/OUTMASK membership pass."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_math import (
    toads_inverse_shape_prior_energy,
    toads_mask_intensity_energy,
    toads_object_intensity_energy,
    toads_outlier_intensity_energy,
    toads_outlier_membership_mix,
    toads_shape_relation_factors,
)


@dataclass(frozen=True)
class ToadsOutmaskMembershipPass:
    """Data maps, updated class/auxiliary memberships and objective."""

    class_data_energy: NDArray[np.float32]
    outlier_data_energy: NDArray[np.float32]
    memberships: NDArray[np.float32]
    outlier_membership: NDArray[np.float32]
    objective: np.float32


def toads_outmask_membership_pass(
    image_values: ArrayLike,
    initial_memberships: ArrayLike,
    initial_outlier_membership: ArrayLike,
    class_types: Sequence[str],
    centroids: ArrayLike,
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
) -> ToadsOutmaskMembershipPass:
    """Compose fixed-input OBJECT/MASK/one-OUTMASK data and update terms.

    The OUTMASK data term is source `maskDistance(image, 0)` and does not use
    the optional field or PSF bias. Its auxiliary term is the weighted fixed
    outlier scalar. The selected source's unbraced LESION smoothing chain
    skips smoothing on that auxiliary energy, while the ordinary OUTMASK data
    term gets smoothing. The separately braced prior path adds its prior to
    both energies. A local sweep preserves those branches instead of reusing
    the otherwise similar OUTLIER sweep, whose auxiliary smoothing differs.

    All other maps and statistics are explicit, supplied inputs in
    `(x,y,z,modality/class)` order. All classes must have active, nonzero
    interior shapes; the source leaves suppressed-class objective energy
    stale, which is outside this focused fixed-input composition.
    """

    types = tuple(str(kind).upper() for kind in class_types)
    if (
        not types
        or types.count("OUTMASK") != 1
        or not set(types).intersection({"OBJECT", "MASK"})
        or any(kind not in {"OBJECT", "MASK", "OUTMASK"} for kind in types)
    ):
        raise ValueError("class_types must contain OBJECT or MASK and exactly one OUTMASK class")

    images = np.asarray(image_values, dtype=np.float32)
    initial = np.asarray(initial_memberships, dtype=np.float32)
    initial_outlier = np.asarray(initial_outlier_membership, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    labels = np.asarray(closest_labels)
    if images.ndim != 4 or images.shape[-1] < 1:
        raise ValueError("image_values must have shape (x,y,z,modality)")
    spatial = images.shape[:3]
    channels, classes = images.shape[-1], len(types)
    if initial.shape != spatial + (classes,):
        raise ValueError("initial_memberships must match the class grid")
    if initial_outlier.shape != spatial:
        raise ValueError("initial_outlier_membership must match the spatial grid")
    if centers.shape != (channels, classes):
        raise ValueError("centroids must match (modality,class) dimensions")
    if shapes.shape != initial.shape:
        raise ValueError("transformed_shapes must match initial_memberships")
    if labels.shape != spatial + (4,):
        raise ValueError("closest_labels must match the spatial grid with four slots")
    if np.any(shapes[1:-1, 1:-1, 1:-1] == np.float32(0.0)):
        raise ValueError("all classes must have active nonzero interior shapes")

    object_energy = toads_object_intensity_energy(
        images, centers, modality_weights, class_variance_weights, image_variances,
        intensity_field=intensity_field, psf_bias=psf_bias,
    )
    mask_energy = toads_mask_intensity_energy(
        images, modality_weights, class_variance_weights, image_variances
    )
    _, outlier_by_class = toads_outlier_intensity_energy(
        images, centers, modality_weights, class_variance_weights, image_variances,
        outlier_value, intensity_field=intensity_field, psf_bias=psf_bias,
    )
    class_energy = np.zeros(spatial + (classes,), dtype=np.float32)
    auxiliary_energy = np.zeros(spatial, dtype=np.float32)
    interior = (slice(1, -1), slice(1, -1), slice(1, -1))
    outmask_index = types.index("OUTMASK")
    for class_index, kind in enumerate(types):
        if kind == "OBJECT":
            class_energy[interior + (class_index,)] = object_energy[interior + (class_index,)]
        else:
            class_energy[interior + (class_index,)] = mask_energy[interior + (class_index,)]
    auxiliary_energy[interior] = outlier_by_class[interior + (outmask_index,)]

    prior_inverse, prior_objective = toads_inverse_shape_prior_energy(
        shapes, atlas_shape_priors, centroid_similarity, prior_coefficient
    )
    factors, suppressed = toads_shape_relation_factors(shapes, labels, prior_factor)
    if np.any(suppressed[1:-1, 1:-1, 1:-1]):
        raise ValueError("all classes must have active nonzero interior shapes")
    if not np.isfinite(class_energy).all() or not np.isfinite(auxiliary_energy).all():
        raise ValueError("OUTMASK class and auxiliary energies must be finite")
    if not np.isfinite(initial).all() or not np.isfinite(initial_outlier).all():
        raise ValueError("initial memberships must be finite")

    updated = initial.copy()
    outlier = initial_outlier.copy()
    smooth = np.float32(smoothing)
    zero = np.float32(1e-10)
    inf = np.float32(1e10)
    if not np.isfinite(smooth) or smooth < 0:
        raise ValueError("smoothing must be finite and nonnegative")
    total_objective = np.float32(0.0)
    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    nx, ny, nz = spatial
    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                denominator = np.float32(0.0)
                objective_at_voxel = np.zeros(classes, dtype=np.float32)
                for class_index, kind in enumerate(types):
                    num = class_energy[x, y, z, class_index]
                    num_outlier = auxiliary_energy[x, y, z]
                    objective_value = num
                    current_outlier = outlier[x, y, z]
                    if kind == "OUTMASK":
                        objective_value = np.float32(
                            np.float32(np.float32(1.0) - current_outlier) * num
                            + np.float32(current_outlier * num_outlier)
                        )

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
                        smoothing_increment = np.float32(
                            np.float32(smoothing_priors[class_index]) * smooth
                        )
                        smoothing_increment = np.float32(smoothing_increment * neighbor_sum)
                        num = np.float32(num + smoothing_increment)
                        # In the selected Java source the OUTMASK `numO`
                        # smoothing addition is unreachable through its
                        # dangling-else control flow. Keep the objective's
                        # generic half-smoothing term regardless.
                        smoothing_objective = np.float32(
                            np.float32(np.float32(0.5) * np.float32(smoothing_priors[class_index]))
                            * smooth
                        )
                        smoothing_objective = np.float32(smoothing_objective * neighbor_sum)
                        objective_value = np.float32(objective_value + smoothing_objective)

                    if np.float32(prior_coefficient) > np.float32(0.0):
                        num = np.float32(num + prior_inverse[x, y, z, class_index])
                        if kind == "OUTMASK":
                            num_outlier = np.float32(
                                num_outlier + prior_inverse[x, y, z, class_index]
                            )
                        objective_value = np.float32(
                            objective_value + prior_objective[x, y, z, class_index]
                        )

                    class_reciprocal = np.float32(1.0) / num if num > zero else inf
                    if kind == "OUTMASK":
                        outlier_reciprocal = (
                            np.float32(1.0) / num_outlier if num_outlier > zero else inf
                        )
                        value, outlier[x, y, z] = toads_outlier_membership_mix(
                            class_reciprocal,
                            current_outlier,
                            outlier_reciprocal,
                            factors[x, y, z, class_index],
                        )
                    else:
                        value = np.float32(factors[x, y, z, class_index] * class_reciprocal)
                    updated[x, y, z, class_index] = value
                    denominator = np.float32(denominator + value)
                    objective_at_voxel[class_index] = objective_value

                if denominator > np.float32(0.0):
                    for class_index in range(classes):
                        updated[x, y, z, class_index] = np.float32(
                            updated[x, y, z, class_index] / denominator
                        )
                        if class_index == outmask_index:
                            outlier[x, y, z] = np.float32(outlier[x, y, z] / denominator)
                        contribution = np.float32(
                            objective_at_voxel[class_index] * updated[x, y, z, class_index]
                        )
                        contribution = np.float32(
                            contribution * updated[x, y, z, class_index]
                        )
                        total_objective = np.float32(total_objective + contribution)

    return ToadsOutmaskMembershipPass(
        class_energy, auxiliary_energy, updated, outlier, total_objective
    )
