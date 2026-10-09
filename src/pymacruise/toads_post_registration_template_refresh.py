"""Compose the selected TOADS post-registration recrop/template refresh.

This is the bounded wrapper sequence after initial rigid registration in
``AlgorithmMultipleToads.runAlgorithm``: uncrop the image, recompute bounds,
crop with the registered transform, then resample atlas labels onto the new
crop. It does not implement or validate TOADS segmentation as a whole.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``, describe
atlas registration and topology priors as parts of TOADS. This adapter
rebuilds the source wrapper's cropped image/template state after rigid
registration so later stages can use the registered prior. The uncrop,
recropping, transform restoration, and resampling sequence is selected
VandyCRUISE wrapper behavior, not a paper-specified algorithm or a complete
registration/segmentation validation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_crop import (
    crop_toads_images,
    find_toads_crop_bounds,
    transform_toads_topology,
    uncrop_toads_images,
)


@dataclass(frozen=True)
class ToadsPostRegistrationTemplateRefresh:
    """State produced by the selected post-registration wrapper sequence."""

    native_images: NDArray[np.float32]
    recropped_images: NDArray[np.float32]
    registered_transform: NDArray[np.float32]
    restored_transform: NDArray[np.float32]
    recropped_transform: NDArray[np.float32]
    bounds: tuple[int, int, int, int, int, int]
    transformed_topology: NDArray[Any]
    call_order: tuple[str, ...]


POST_REGISTRATION_REFRESH_CALL_ORDER = (
    "uncropImages",
    "findCroppingBoundaries",
    "cropImages",
    "setImageInfo",
    "computeTransformedTemplate",
    "updateTransformedTemplate",
)


def toads_post_registration_template_refresh(
    *,
    cropped_images: Any,
    native_shape_xyz: tuple[int, int, int],
    previous_bounds: tuple[int, int, int, int, int, int],
    padding: int,
    registered_crop_transform: Any,
    topology_template: Any,
    image_resolution_mm: tuple[float, float, float],
    template_resolution_mm: tuple[float, float, float],
    orientation_signs: tuple[float, float, float],
    robust_ranges: Any,
    background_ratio: float,
    background_label: int = 1,
) -> ToadsPostRegistrationTemplateRefresh:
    """Run the crop/template refresh following rigid atlas registration.

    Arrays use the selected ``[channel,x,y,z]`` / ``[x,y,z]`` convention. The
    supplied ``registered_crop_transform`` is the shared preprocessing
    transform after ``registerShapesPyramid`` mutated it. Bounds are recomputed
    on the restored native image grid using the still-template-space topology;
    the resulting labels are then sampled over the *entire* recropped grid,
    including the padding voxels, as ``MultipleToadDeformableAtlas`` does.
    Input arrays are treated as read-only.
    """
    images = np.asarray(cropped_images, dtype=np.float32)
    template = np.asarray(topology_template)
    transform = np.asarray(registered_crop_transform, dtype=np.float32)
    ranges = np.asarray(robust_ranges, dtype=np.float32)
    if images.ndim != 4 or any(size < 1 for size in images.shape):
        raise ValueError("cropped_images must have shape (channels,x,y,z)")
    if template.ndim != 3 or any(size < 2 for size in template.shape):
        raise ValueError("topology_template must be a 3D template label array")
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("registered_crop_transform must contain six finite values")
    if ranges.shape != (images.shape[0], 2):
        raise ValueError("robust_ranges must have shape (channels,2)")
    if len(native_shape_xyz) != 3 or min(native_shape_xyz) < 1:
        raise ValueError("native_shape_xyz must contain three positive dimensions")
    if padding < 0:
        raise ValueError("padding must be nonnegative")

    native_images, restored_transform = uncrop_toads_images(
        images,
        native_shape_xyz,
        previous_bounds,
        padding,
        transform,
        image_resolution=image_resolution_mm,
        orientation_signs=orientation_signs,
    )
    bounds = find_toads_crop_bounds(
        native_images,
        template,
        restored_transform,
        image_resolution=image_resolution_mm,
        template_resolution=template_resolution_mm,
        orientation_signs=orientation_signs,
        robust_ranges=ranges,
        background_ratio=background_ratio,
        background_label=background_label,
    )
    cropped = crop_toads_images(
        native_images,
        bounds,
        padding,
        restored_transform,
        image_resolution=image_resolution_mm,
        orientation_signs=orientation_signs,
    )
    full_crop_bounds = (
        0,
        cropped.images.shape[1] - 1,
        0,
        cropped.images.shape[2] - 1,
        0,
        cropped.images.shape[3] - 1,
    )
    transformed = transform_toads_topology(
        template,
        full_crop_bounds,
        0,
        cropped.transform,
        crop_shape=tuple(int(size) for size in cropped.images.shape[1:]),
        image_resolution=image_resolution_mm,
        template_resolution=template_resolution_mm,
        orientation_signs=orientation_signs,
        # AlgorithmMultipleToads changes atlas image info after crop, but the
        # ParametricTransform created before registration retains the original
        # cropped dimensions used for its x0i/y0i/z0i center fields.
        image_center_shape=tuple(int(size) for size in images.shape[1:]),
        background_label=background_label,
    )
    return ToadsPostRegistrationTemplateRefresh(
        native_images=native_images,
        recropped_images=cropped.images,
        registered_transform=transform.copy(),
        restored_transform=restored_transform,
        recropped_transform=cropped.transform,
        bounds=bounds,
        transformed_topology=transformed,
        call_order=POST_REGISTRATION_REFRESH_CALL_ORDER,
    )
