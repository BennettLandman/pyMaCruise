"""Fixed-input exact-membership dispatch for ordinary and optimized TOADS types.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``. The
source-selected optimized tissue energy is one statistical classification
component, not a complete paper-level segmentation implementation.
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
    toads_ordered_membership_sweep,
)


_STANDARD_TYPES = {"OBJECT", "MASK"}
_OPTIMIZED_TYPES = {"OPTIMIZED_CSF", "OPTIMIZED_GM", "OPTIMIZED_WM"}


@dataclass(frozen=True)
class ToadsOptimizedDispatchResult:
    """Updated class memberships and objective for one fixed-input pass."""

    memberships: NDArray[np.float32]
    objective: np.float32


def toads_optimized_membership_dispatch(
    image_values: ArrayLike,
    initial_memberships: ArrayLike,
    class_types: Sequence[str],
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
    intensity_field: ArrayLike | None = None,
    psf_bias: ArrayLike | None = None,
) -> ToadsOptimizedDispatchResult:
    """Compose OBJECT/MASK and optimized tissue intensity branches.

    Supported type lists contain any ordered mixture of OBJECT, MASK,
    OPTIMIZED_CSF, OPTIMIZED_GM, and OPTIMIZED_WM, and must include at least
    one optimized tissue class. Centroids, each optimized class's cached
    PowerTable exponent, variance weights, normalized image-variance weights,
    shape samples, similarity, and priors are explicit fixed inputs. This
    function computes no statistics and composes only the source intensity
    terms with the validated shared ordered membership sweep.

    Optimized likelihoods use the selected class-specific exponent matrix;
    OBJECT/MASK energy maps are selected only for their corresponding class
    columns. The existing sweep preserves class ordering, in-place
    six-neighbor smoothing, reciprocal/relation update, normalization,
    objective accumulation, and untouched borders.

    Inputs are already aligned in RAS voxel order and use ``(x,y,z,channel)``
    for image-like arrays and ``(x,y,z,class)`` for class-like arrays. LESION,
    OUTLIER, OUTMASK, RAYLEIGH, auxiliary maps, and mixtures involving them
    are rejected because their special data and membership-map branches are
    not part of this composition.
    """

    types = tuple(str(value).upper() for value in class_types)
    if not types:
        raise ValueError("class_types must contain at least one class")
    unknown = sorted(set(types) - (_STANDARD_TYPES | _OPTIMIZED_TYPES))
    if unknown:
        raise ValueError(f"unsupported optimized-dispatch class types: {', '.join(unknown)}")
    if not set(types) & _OPTIMIZED_TYPES:
        raise ValueError("class_types must include at least one optimized tissue class")

    initial = np.asarray(initial_memberships)
    if initial.ndim != 4 or len(types) != initial.shape[-1]:
        raise ValueError("class_types must name every membership class in order")

    object_energy = toads_object_intensity_energy(
        image_values,
        centroids,
        modality_weights,
        class_variance_weights,
        image_variances,
        intensity_field=intensity_field,
        psf_bias=psf_bias,
    )
    mask_energy = toads_mask_intensity_energy(
        image_values, modality_weights, class_variance_weights, image_variances
    )
    energies = np.empty_like(object_energy)
    centers = np.asarray(centroids, dtype=np.float32)
    exponents = np.asarray(exponent_by_channel_class, dtype=np.float32)
    variance_weights = np.asarray(class_variance_weights, dtype=np.float32)
    if exponents.shape != centers.shape:
        raise ValueError("exponent_by_channel_class must match centroid dimensions")

    for index, kind in enumerate(types):
        if kind == "OBJECT":
            energies[..., index] = object_energy[..., index]
        elif kind == "MASK":
            energies[..., index] = mask_energy[..., index]
        else:
            # Run the validated PowerTable likelihood on this category alone,
            # so unused OBJECT/MASK exponent columns cannot affect its branch.
            class_energy = toads_optimized_tissue_data_term(
                image_values,
                centers[:, index : index + 1],
                exponents[:, index : index + 1],
                modality_weights,
                variance_weights[:, index : index + 1],
                image_variances,
                intensity_field=intensity_field,
                psf_bias=psf_bias,
            )
            energies[..., index] = class_energy[..., 0]

    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    inverse_prior, prior_objective = toads_inverse_shape_prior_energy(
        shapes, atlas_shape_priors, centroid_similarity, prior_coefficient
    )
    updated, objective = toads_ordered_membership_sweep(
        energies,
        initial_memberships,
        shapes,
        closest_labels,
        smoothing_priors,
        smoothing,
        prior_factor,
        prior_inverse_energy=inverse_prior,
        prior_objective_energy=prior_objective,
    )
    return ToadsOptimizedDispatchResult(updated, objective)
