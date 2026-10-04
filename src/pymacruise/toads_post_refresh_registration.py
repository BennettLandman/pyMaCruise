"""Compose exact TOADS memberships and the second rigid registration update.

The selected wrapper reinitializes memberships after its post-pyramid crop and
then performs a single-scale `registerShapes()` update. This adapter composes
those calls from the existing source-matched initializer and registration
kernels; it is not TOADS segmentation iteration.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_initialization import initial_exact_memberships
from .toads_registration_sampling import (
    ToadsRigidRegistrationResult,
    register_toads_atlas_rigid,
)
from .toads_shapes import ToadShapeBounds, _rotation_matrix, compute_toad_shape_bounds


@dataclass(frozen=True)
class ToadsPostRefreshRegistration:
    """Reinitialized exact memberships and the selected single-scale update."""

    memberships: NDArray[np.float32]
    registration: ToadsRigidRegistrationResult
    final_transform: NDArray[np.float32]
    final_rotation: NDArray[np.float32]
    final_shape_bounds: tuple[ToadShapeBounds, ...]
    registered_classes: NDArray[np.bool_]
    call_order: tuple[str, ...]


def register_toads_after_template_refresh(
    *,
    recropped_images: Any,
    centroids: Any,
    atlas_shapes: Any,
    initial_transform: Any,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    background_class_index: int,
    crop_border: int = 10,
    max_iterations: int = 5,
    energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
) -> ToadsPostRefreshRegistration:
    """Reinitialize exact memberships, then run the selected rigid update.

    `recropped_images` must be channel-first and already in RAS voxel order.
    Atlas maps remain in template-grid RAS order. The source wrapper has
    `useApprox=false`, so this deliberately selects exact intensity-only
    initialization; the optional atlas argument in the Java overload is not
    used by that exact branch. Registration is the selected scale-one
    `registerShapes()` call with subsample 2 and source transformed-shape
    bounds. The caller controls iteration count for focused synthetic runs.
    """
    images = np.asarray(recropped_images, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    transform = np.asarray(initial_transform, dtype=np.float32)
    if images.ndim != 4 or any(size <= 0 for size in images.shape):
        raise ValueError("recropped_images must have shape (channels,x,y,z)")
    if centers.ndim != 2 or centers.shape[0] != images.shape[0]:
        raise ValueError("centroids must have shape (channels,classes)")
    if shapes.ndim != 4 or shapes.shape[-1] != centers.shape[1]:
        raise ValueError("atlas_shapes must be a template grid with one map per centroid class")
    if not np.isfinite(images).all() or not np.isfinite(centers).all() or not np.isfinite(shapes).all():
        raise ValueError("images, centroids, and atlas shapes must be finite")
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_transform must contain six finite rigid parameters")
    if not isinstance(background_class_index, (int, np.integer)) or not 0 <= background_class_index < centers.shape[1]:
        raise ValueError("background_class_index must identify one atlas class")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, (int, np.integer)) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")

    memberships = initial_exact_memberships(images, centers, crop_border=crop_border)
    bounds = tuple(
        compute_toad_shape_bounds(
            shapes[..., klass],
            image_dimensions_xyz=tuple(int(v) for v in memberships.shape[:3]),
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=tuple(float(v) for v in transform),
        )
        for klass in range(centers.shape[1])
    )
    low = np.asarray([item.lower_xyz for item in bounds], dtype=np.int32)
    high = np.asarray([item.upper_xyz for item in bounds], dtype=np.int32)
    registered = np.ones(centers.shape[1], dtype=bool)
    registered[background_class_index] = False
    registration = register_toads_atlas_rigid(
        memberships,
        shapes,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        initial_transform=tuple(float(v) for v in transform),
        class_bounds_low_xyz=low,
        class_bounds_high_xyz=high,
        registered_classes=registered,
        offset=0,
        subsample=2,
        max_iterations=int(max_iterations),
        energy_tolerance=energy_tolerance,
        minimum_lambda=minimum_lambda,
        # MultipleToadDeformableAtlas defaults precompute=true in the selected
        # VandyCRUISE source; use the same affine coordinate path here.
        precomputed_affine=True,
    )
    final_bounds = tuple(
        compute_toad_shape_bounds(
            shapes[..., klass],
            image_dimensions_xyz=tuple(int(v) for v in memberships.shape[:3]),
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=tuple(float(v) for v in registration.transform),
        )
        for klass in range(centers.shape[1])
    )
    return ToadsPostRefreshRegistration(
        memberships,
        registration,
        registration.transform.copy(),
        _rotation_matrix(registration.transform[:3]),
        final_bounds,
        registered,
        (
            "segmentation.initCentroids",
            "preprocess.initialMemberships",
            "atlas.updateShapeRegistration",
            "atlas.registerShapes()",
        ),
    )
