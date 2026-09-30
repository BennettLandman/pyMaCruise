"""Selected base-TOADS exact-membership dispatch over validated kernels.

This module is deliberately narrower than
``MultipleToadSegmentation.computeExactMemberships``: it dispatches the
source-supported ordinary OBJECT/MASK[/OUTLIER] branches into their existing
fixed-input kernels. Image statistics, class parameters, priors, relation
labels, and auxiliary-map initialization are explicit inputs. It does not
estimate statistics or claim to implement the complete TOADS loop. Paper
context: Bazin and Pham (2007, 2008), see ``docs/papers/toads.md``; this is
only a source-selected statistical-classification substep.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_math import (
    toads_object_mask_membership_pass,
    toads_object_mask_outlier_membership_pass,
)


@dataclass(frozen=True)
class ToadsExactMembershipDispatch:
    """Result of one selected exact-membership dispatch pass.

    ``outlier_membership`` is present only for a dispatched OUTLIER class.
    The objective is the selected source's float32 accumulated objective.
    """

    memberships: NDArray[np.float32]
    objective: np.float32
    outlier_membership: NDArray[np.float32] | None = None


def toads_exact_membership_dispatch(
    image_values: ArrayLike,
    initial_memberships: ArrayLike,
    class_types: Sequence[str],
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
    initial_outlier_membership: ArrayLike | None = None,
    outlier_value: float | None = None,
    intensity_field: ArrayLike | None = None,
    psf_bias: ArrayLike | None = None,
) -> ToadsExactMembershipDispatch:
    """Dispatch a fixed-statistics pass for supported base TOADS class mixes.

    Supported ordered type lists are OBJECT/MASK mixtures, or such a mixture
    plus exactly one OUTLIER. Class order is preserved. The existing helpers
    provide the source-matched data terms, prior construction, asynchronous
    six-face update, auxiliary outlier-map update, normalization, and
    objective accumulation.

    All spatial arrays use ``(x, y, z[, modality/class])``; no reorientation is
    performed here because these ndarray kernels consume already aligned
    grids. The caller must supply inputs in RAS voxel order. The only supported
    auxiliary map is the scalar OUTLIER membership map.

    LESION, OUTMASK, RAYLEIGH, optimized tissue types, and any mixtures that
    contain those types are rejected. Their source branches require different
    class-specific and auxiliary-map semantics and are not silently treated
    as OBJECT. The source atlas parser's unknown-name-to-OBJECT fallback is
    likewise outside this typed Python API.
    """

    types = tuple(str(value).upper() for value in class_types)
    if not types:
        raise ValueError("class_types must contain at least one class")
    supported = {"OBJECT", "MASK", "OUTLIER"}
    unknown = sorted(set(types) - supported)
    if unknown:
        raise ValueError(f"unsupported exact-membership class types: {', '.join(unknown)}")
    if types.count("OUTLIER") > 1:
        raise ValueError("class_types may contain at most one OUTLIER class")
    if len(types) != np.asarray(initial_memberships).shape[-1]:
        raise ValueError("class_types must name every membership class in order")

    has_outlier = "OUTLIER" in types
    if has_outlier:
        if initial_outlier_membership is None or outlier_value is None:
            raise ValueError("OUTLIER dispatch requires its initial map and outlier_value")
        memberships, outlier, objective = toads_object_mask_outlier_membership_pass(
            image_values,
            initial_memberships,
            initial_outlier_membership,
            centroids,
            modality_weights,
            class_variance_weights,
            image_variances,
            outlier_value,
            transformed_shapes,
            closest_labels,
            smoothing_priors,
            centroid_similarity,
            atlas_shape_priors,
            smoothing,
            prior_factor,
            prior_coefficient,
            class_types=types,
            intensity_field=intensity_field,
            psf_bias=psf_bias,
        )
        return ToadsExactMembershipDispatch(memberships, objective, outlier)

    if initial_outlier_membership is not None or outlier_value is not None:
        raise ValueError("OUTLIER inputs are only valid when class_types contains OUTLIER")
    memberships, objective = toads_object_mask_membership_pass(
        image_values,
        initial_memberships,
        centroids,
        modality_weights,
        class_variance_weights,
        image_variances,
        transformed_shapes,
        closest_labels,
        smoothing_priors,
        centroid_similarity,
        atlas_shape_priors,
        smoothing,
        prior_factor,
        prior_coefficient,
        class_types=types,
        intensity_field=intensity_field,
        psf_bias=psf_bias,
    )
    return ToadsExactMembershipDispatch(memberships, objective)
