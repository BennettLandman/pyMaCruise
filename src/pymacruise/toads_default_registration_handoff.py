"""Compose the selected default TOADS membership-to-rigid-registration handoff."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_default_atlas_membership import (
    ToadsDefaultAtlasMembershipPass,
    toads_default_atlas_membership_pass,
)
from .toads_registration_sampling import (
    ToadsRigidRegistrationResult,
    register_toads_atlas_rigid,
)
from .toads_shapes import ToadShapeBounds, _rotation_matrix, compute_toad_shape_bounds


@dataclass(frozen=True)
class ToadsDefaultRegistrationHandoff:
    """Result after default membership update and optional rigid registration."""

    membership_pass: ToadsDefaultAtlasMembershipPass
    registration: ToadsRigidRegistrationResult | None
    transform: NDArray[np.float32]
    rotation: NDArray[np.float32]
    shape_bounds: tuple[ToadShapeBounds, ...]


def _bounds_for_transform(
    atlas_shapes: NDArray[np.float32],
    image_dimensions_xyz: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: NDArray[np.float32],
) -> tuple[ToadShapeBounds, ...]:
    return tuple(
        compute_toad_shape_bounds(
            atlas_shapes[..., class_index],
            image_dimensions_xyz=image_dimensions_xyz,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=tuple(float(value) for value in transform),
        )
        for class_index in range(atlas_shapes.shape[-1])
    )


def toads_default_atlas_registration_handoff(
    *,
    image_values: Any,
    initial_memberships: Any,
    centroids: Any,
    exponent_by_channel_class: Any,
    modality_weights: Any,
    class_variance_weights: Any,
    image_variances: Any,
    transformed_shapes: Any,
    closest_labels: Any,
    smoothing_priors: Any,
    centroid_similarity: Any,
    atlas_shape_priors: Any,
    smoothing: float,
    prior_factor: float,
    prior_coefficient: float,
    atlas_shapes: Any,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    initial_transform: tuple[float, float, float, float, float, float],
    align_shapes: bool = True,
    intensity_field: Any | None = None,
    psf_bias: Any | None = None,
    offset: int = 0,
    subsample: int = 1,
    max_iterations: int = 50,
    energy_tolerance: float = 1e-4,
    minimum_lambda: float = 1e-3,
) -> ToadsDefaultRegistrationHandoff:
    """Run one default 11-class pass, then the wrapper's optional `registerShapes`.

    The `transformed_shapes` input is the current-transform prior volume used
    by `computeExactMemberships`. `atlas_shapes` is the template-space atlas
    passed to the rigid registration optimizer. Registration consumes the
    just-updated membership maps. The wrapper invokes this registration before
    its convergence/iteration-limit checks. Stop-path growth and continuing
    relation refresh are outside this composition.
    """
    membership = toads_default_atlas_membership_pass(
        image_values=image_values,
        initial_memberships=initial_memberships,
        centroids=centroids,
        exponent_by_channel_class=exponent_by_channel_class,
        modality_weights=modality_weights,
        class_variance_weights=class_variance_weights,
        image_variances=image_variances,
        transformed_shapes=transformed_shapes,
        closest_labels=closest_labels,
        smoothing_priors=smoothing_priors,
        centroid_similarity=centroid_similarity,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=smoothing,
        prior_factor=prior_factor,
        prior_coefficient=prior_coefficient,
        intensity_field=intensity_field,
        psf_bias=psf_bias,
    )
    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    params = np.asarray(initial_transform, dtype=np.float32)
    if shapes.ndim != 4 or shapes.shape[-1] != 11 or min(shapes.shape[:3]) < 2:
        raise ValueError("atlas_shapes must be a 4D template grid with eleven classes")
    if not np.isfinite(shapes).all():
        raise ValueError("atlas_shapes must be finite")
    if params.shape != (6,) or not np.isfinite(params).all():
        raise ValueError("initial_transform must contain six finite rigid parameters")

    current_bounds = _bounds_for_transform(
        shapes, membership.memberships.shape[:3], image_resolutions_mm,
        template_resolutions_mm, params
    )
    if not align_shapes:
        return ToadsDefaultRegistrationHandoff(
            membership,
            None,
            params.copy(),
            _rotation_matrix(params[:3]),
            current_bounds,
        )

    low = np.asarray([item.lower_xyz for item in current_bounds], dtype=np.int32)
    high = np.asarray([item.upper_xyz for item in current_bounds], dtype=np.int32)
    registered = np.ones(11, dtype=bool)
    registered[0] = False
    registration = register_toads_atlas_rigid(
        membership.memberships,
        shapes,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        initial_transform=tuple(float(value) for value in params),
        class_bounds_low_xyz=low,
        class_bounds_high_xyz=high,
        registered_classes=registered,
        offset=offset,
        subsample=subsample,
        max_iterations=max_iterations,
        energy_tolerance=energy_tolerance,
        minimum_lambda=minimum_lambda,
    )
    final_bounds = _bounds_for_transform(
        shapes, membership.memberships.shape[:3], image_resolutions_mm,
        template_resolutions_mm, registration.transform
    )
    return ToadsDefaultRegistrationHandoff(
        membership,
        registration,
        registration.transform.copy(),
        _rotation_matrix(registration.transform[:3]),
        final_bounds,
    )
