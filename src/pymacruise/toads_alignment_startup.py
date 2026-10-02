"""Compose the selected pre-registration TOADS alignment startup stages.

This is deliberately a bounded adapter around already source-checked pieces:
mask/normalize the native T1, find/crop image and topology support, initialize
intensity centroids, and compute the exact-membership array supplied to rigid
atlas registration. It does not perform atlas registration or TOADS evolution.
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib
import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_crop import crop_toads_images, find_toads_crop_bounds, transform_toads_topology
from .toads_initialization import initial_exact_memberships, initial_intensity_centroids
from .toads_preprocessing import ToadsIntensityResult, ToadsRASGeometry, prepare_toads_intensity, toads_ras_geometry


@dataclass(frozen=True)
class ToadsAlignmentStartup:
    """Source-ordered state immediately before the first rigid atlas update."""

    intensity: ToadsIntensityResult
    geometry: ToadsRASGeometry
    crop_bounds: tuple[int, int, int, int, int, int]
    cropped_images: NDArray[np.float32]
    cropped_topology: NDArray
    transform: NDArray[np.float32]
    centroids: NDArray[np.float32]
    initial_memberships: NDArray[np.float32]
    stage_order: tuple[str, ...]


def prepare_toads_alignment_startup(
    t1_n4: nib.spatialimages.SpatialImage,
    brain_mask: nib.spatialimages.SpatialImage,
    topology_ras: ArrayLike,
    intensity_priors: ArrayLike,
    template_resolution_mm: tuple[float, float, float],
    *,
    padding: int = 10,
    background_label: int = 1,
    background_ratio: float = 0.0,
) -> ToadsAlignmentStartup:
    """Prepare one-channel base TOADS state up to rigid atlas registration.

    ``topology_ras`` is an atlas label array already loaded/reoriented into the
    selected RAS voxel-axis order; this function cannot infer orientation from
    an array without an affine. ``template_resolution_mm`` comes from the
    bundled atlas descriptor, which the selected Java loader uses in place of
    any conflicting NIfTI header spacing. The T1 and user mask are normalized
    to RAS by :func:`prepare_toads_intensity` and non-RAS inputs warn there.

    The selected wrapper initializes atlas alignment from exact memberships.
    This adapter stops at that registration input and does not claim an atlas
    fit, complete TOADS iteration, or segmentation equivalence.
    """
    intensity = prepare_toads_intensity(t1_n4, brain_mask)
    geometry = toads_ras_geometry(intensity.image)
    values = np.asarray(geometry.image.dataobj, dtype=np.float32)
    topology = np.asarray(topology_ras)
    priors = np.asarray(intensity_priors, dtype=np.float32)
    template_resolution = np.asarray(template_resolution_mm, dtype=np.float32)

    if values.ndim != 3 or topology.ndim != 3 or topology.size == 0:
        raise ValueError("T1 and RAS atlas topology must be nonempty 3D arrays")
    if not np.issubdtype(topology.dtype, np.integer):
        raise ValueError("RAS atlas topology must contain integer structure labels")
    if priors.ndim != 2 or priors.shape[0] != 1 or priors.shape[1] < 1:
        raise ValueError("base TOADS intensity priors must have shape (1, classes)")
    if not np.isfinite(priors).all():
        raise ValueError("base TOADS intensity priors must be finite")
    if template_resolution.shape != (3,) or not np.isfinite(template_resolution).all() or np.any(template_resolution <= 0):
        raise ValueError("template_resolution_mm must contain three finite positive values")
    if not isinstance(padding, (int, np.integer)) or padding < 0:
        raise ValueError("padding must be a non-negative integer")

    image_channels = values[np.newaxis, ...]
    initial_transform = np.asarray(geometry.initial_transform, dtype=np.float32)
    robust_ranges = np.asarray(
        [[intensity.robust_minimum, intensity.robust_maximum]], dtype=np.float32
    )
    bounds = find_toads_crop_bounds(
        image_channels,
        topology,
        initial_transform,
        image_resolution=tuple(float(v) for v in np.abs(geometry.signed_image_resolutions_mm)),
        template_resolution=tuple(float(v) for v in template_resolution),
        orientation_signs=geometry.orientation_signs,
        robust_ranges=robust_ranges,
        background_ratio=background_ratio,
        background_label=background_label,
        normalized=True,
    )
    crop = crop_toads_images(
        image_channels,
        bounds,
        int(padding),
        initial_transform,
        image_resolution=tuple(float(v) for v in np.abs(geometry.signed_image_resolutions_mm)),
        orientation_signs=geometry.orientation_signs,
    )
    cropped_topology = transform_toads_topology(
        topology,
        bounds,
        int(padding),
        crop.transform,
        crop_shape=tuple(int(v) for v in crop.images.shape[1:]),
        image_resolution=tuple(float(v) for v in np.abs(geometry.signed_image_resolutions_mm)),
        template_resolution=tuple(float(v) for v in template_resolution),
        orientation_signs=geometry.orientation_signs,
        background_label=background_label,
    )
    centroids = initial_intensity_centroids(
        [intensity.robust_minimum], [intensity.robust_maximum], priors
    )
    memberships = initial_exact_memberships(crop.images, centroids, crop_border=int(padding))
    return ToadsAlignmentStartup(
        intensity=intensity,
        geometry=geometry,
        crop_bounds=bounds,
        cropped_images=crop.images,
        cropped_topology=cropped_topology,
        transform=crop.transform,
        centroids=centroids,
        initial_memberships=memberships,
        stage_order=(
            "normalizeImages",
            "findCroppingBoundaries",
            "cropImages",
            "initialCentroids",
            "initialMemberships",
        ),
    )
