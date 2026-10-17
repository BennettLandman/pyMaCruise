"""Selected input normalization and FCM initialization for AtlasEM."""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Sequence
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.atlas_em_statistics import normalize_atlas_em_priors
from pymacruise.mipav_morphology import mipav_fill_holes_25d
from pymacruise.toads_preprocessing import robust_intensity_bounds
from pymacruise.toads_shape_volume import transform_toads_shape_priors
from pymacruise.toads_shapes import compute_toad_shape_bounds
from pymacruise.toads_shapes import _rotation_matrix


@dataclass(frozen=True)
class AtlasEMNormalizedImages:
    images: NDArray[np.float32]
    robust_minima: NDArray[np.float32]
    robust_maxima: NDArray[np.float32]


@dataclass(frozen=True)
class AtlasEMCroppedImages:
    """Source-shaped AtlasEM crop with inclusive source bounds and transform state."""

    images: NDArray[np.float32]
    lower_xyz: tuple[int, int, int]
    upper_xyz: tuple[int, int, int]
    border: int
    original_dimensions_xyz: tuple[int, int, int]
    signed_resolutions_mm: tuple[float, float, float]
    transform: NDArray[np.float32]


@dataclass(frozen=True)
class AtlasEMDescriptorInitialization:
    """Descriptor-derived initial class means and source inverse-distance memberships."""

    class_names: tuple[str, ...]
    centroids: NDArray[np.float32]
    memberships: NDArray[np.float32]


def build_atlas_em_internal_mask(images: ArrayLike) -> NDArray[np.bool_]:
    """Construct AtlasEM's internal statistical mask from positive channels.

    ``AtlasEMSegmentation`` marks a voxel active only when every input channel
    is strictly positive, then runs MIPAV ``AlgorithmMorphology25D.FILL_HOLES``
    on the resulting mask. This is the selected algorithm's internal mask; it
    does not replace or modify the user-supplied brain mask used by the outer
    MaCRUISE workflow.
    """
    values = np.asarray(images, dtype=np.float32)
    if values.ndim != 4 or any(size < 2 for size in values.shape[:3]) or values.shape[3] < 1:
        raise ValueError("images must be a 4D (x, y, z, channel) grid with nonempty dimensions")
    if not np.isfinite(values).all():
        raise ValueError("images must contain only finite values")
    positive_all_channels = np.all(values > np.float32(0.0), axis=-1)
    return mipav_fill_holes_25d(positive_all_channels)


def find_atlas_em_crop_bounds(
    images: ArrayLike,
    robust_minima: ArrayLike,
    robust_maxima: ArrayLike,
    *,
    background_ratio: float,
    normalized: bool = False,
) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    """Find inclusive image bounds using AtlasEMPreprocess crop predicates."""
    values = np.asarray(images, dtype=np.float32)
    low = np.asarray(robust_minima, dtype=np.float32)
    high = np.asarray(robust_maxima, dtype=np.float32)
    ratio = np.float32(background_ratio)
    if values.ndim != 4 or any(size < 1 for size in values.shape[:3]) or values.shape[3] < 1:
        raise ValueError("images must have shape (x, y, z, nonzero channels)")
    if low.shape != (values.shape[3],) or high.shape != low.shape:
        raise ValueError("robust minima/maxima must match image channels")
    if not np.isfinite(values).all() or not np.isfinite(low).all() or not np.isfinite(high).all() or not np.isfinite(ratio):
        raise ValueError("images, bounds, and background_ratio must be finite")
    thresholds = np.full(values.shape[3], ratio, dtype=np.float32)
    if not normalized:
        thresholds = np.multiply(ratio, np.subtract(high, low, dtype=np.float32), dtype=np.float32)
    foreground = np.any(values > thresholds.reshape((1, 1, 1, -1)), axis=3)
    if not np.any(foreground):
        raise ValueError("AtlasEM crop predicate found no foreground voxels")
    lower: list[int] = []
    upper: list[int] = []
    for axis in range(3):
        present = np.any(foreground, axis=tuple(i for i in range(3) if i != axis))
        lower.append(int(np.argmax(present)))
        upper.append(int(present.size - 1 - np.argmax(present[::-1])))
    return tuple(lower), tuple(upper)


def crop_atlas_em_images(
    images: ArrayLike,
    lower_xyz: Sequence[int],
    upper_xyz: Sequence[int],
    *,
    border: int,
    signed_resolutions_mm: Sequence[float],
    transform: ArrayLike,
) -> AtlasEMCroppedImages:
    """Crop channel images with AtlasEM's zero halo and rigid-translation update."""
    values = np.asarray(images, dtype=np.float32)
    lower = tuple(int(v) for v in lower_xyz)
    upper = tuple(int(v) for v in upper_xyz)
    spacing = np.asarray(signed_resolutions_mm, dtype=np.float32)
    parameters = np.asarray(transform, dtype=np.float32).copy()
    pad = int(border)
    if values.ndim != 4 or values.shape[3] < 1:
        raise ValueError("images must have shape (x, y, z, nonzero channels)")
    if len(lower) != 3 or len(upper) != 3 or any(lo < 0 or hi < lo or hi >= size for lo, hi, size in zip(lower, upper, values.shape[:3], strict=True)):
        raise ValueError("crop bounds must be inclusive and within the image")
    if pad < 0 or spacing.shape != (3,) or not np.isfinite(spacing).all() or np.any(spacing == 0):
        raise ValueError("border must be nonnegative and signed resolutions finite/nonzero")
    if parameters.shape != (6,) or not np.isfinite(parameters).all() or not np.isfinite(values).all():
        raise ValueError("transform must have six finite values and images must be finite")
    cropped_dims = tuple(upper[a] - lower[a] + 1 + 2 * pad for a in range(3))
    cropped = np.zeros(cropped_dims + (values.shape[3],), dtype=np.float32)
    destination = tuple(slice(pad, pad + upper[a] - lower[a] + 1) for a in range(3))
    source = tuple(slice(lower[a], upper[a] + 1) for a in range(3))
    cropped[destination + (slice(None),)] = values[source + (slice(None),)]

    rotation = _rotation_matrix(parameters[:3], java_float32_norm=True)
    centers_delta = np.asarray(
        [
            np.float32(
                np.float32(cropped_dims[a] / 2.0)
                + np.float32(lower[a])
                - np.float32(pad)
                - np.float32(values.shape[a] / 2.0)
            )
            for a in range(3)
        ], dtype=np.float32,
    )
    physical_delta = np.multiply(centers_delta, spacing, dtype=np.float32)
    for axis in range(3):
        rotated = np.float32(0.0)
        for source_axis in range(3):
            rotated = np.float32(rotated + np.float32(rotation[axis, source_axis] * physical_delta[source_axis]))
        parameters[axis + 3] = np.float32(parameters[axis + 3] + rotated)
    return AtlasEMCroppedImages(
        cropped, lower, upper, pad, tuple(values.shape[:3]),
        tuple(float(value) for value in spacing), parameters,
    )


def uncrop_atlas_em_images(crop: AtlasEMCroppedImages) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Restore cropped data and rigid translation using AtlasEM's zero-fill loops."""
    values = np.asarray(crop.images, dtype=np.float32)
    output = np.zeros(crop.original_dimensions_xyz + (values.shape[3],), dtype=np.float32)
    lower, upper, pad = crop.lower_xyz, crop.upper_xyz, crop.border
    destination: list[slice] = []
    source: list[slice] = []
    for axis in range(3):
        low = max(lower[axis] - pad, 0)
        high = min(upper[axis] + pad, output.shape[axis] - 1)
        destination.append(slice(low, high + 1))
        source.append(slice(low - lower[axis] + pad, high - lower[axis] + pad + 1))
    output[tuple(destination) + (slice(None),)] = values[tuple(source) + (slice(None),)]
    parameters = np.asarray(crop.transform, dtype=np.float32).copy()
    rotation = _rotation_matrix(parameters[:3], java_float32_norm=True)
    cropped_dims = values.shape[:3]
    centers_delta = np.asarray(
        [
            np.float32(
                np.float32(cropped_dims[a] / 2.0)
                + np.float32(lower[a])
                - np.float32(pad)
                - np.float32(crop.original_dimensions_xyz[a] / 2.0)
            )
            for a in range(3)
        ], dtype=np.float32,
    )
    physical_delta = np.multiply(
        centers_delta, np.asarray(crop.signed_resolutions_mm, dtype=np.float32), dtype=np.float32
    )
    for axis in range(3):
        rotated = np.float32(0.0)
        for source_axis in range(3):
            rotated = np.float32(rotated + np.float32(rotation[axis, source_axis] * physical_delta[source_axis]))
        parameters[axis + 3] = np.float32(parameters[axis + 3] - rotated)
    return output, parameters


def recrop_atlas_em_after_registration(
    crop: AtlasEMCroppedImages,
    registered_transform: ArrayLike,
    robust_minima: ArrayLike,
    robust_maxima: ArrayLike,
    *,
    background_ratio: float,
) -> AtlasEMCroppedImages:
    """Reproduce AtlasEM's uncrop, recompute-boundaries, and crop sequence.

    The optimizer's current crop-space transform is first installed into the
    crop state, as the Java atlas and preprocessing object share that array.
    The routine then returns to native-grid coordinates, finds the updated
    robust-span foreground box, and applies the requested zero border again.
    """
    current = np.asarray(registered_transform, dtype=np.float32)
    if current.shape != (6,) or not np.isfinite(current).all():
        raise ValueError("registered_transform must contain six finite values")
    registered_crop = replace(crop, transform=current.copy())
    native_images, native_transform = uncrop_atlas_em_images(registered_crop)
    lower, upper = find_atlas_em_crop_bounds(
        native_images,
        robust_minima,
        robust_maxima,
        background_ratio=background_ratio,
    )
    return crop_atlas_em_images(
        native_images,
        lower,
        upper,
        border=crop.border,
        signed_resolutions_mm=crop.signed_resolutions_mm,
        transform=native_transform,
    )


def transform_atlas_em_rigid_priors(
    atlas_shapes: ArrayLike,
    mask: ArrayLike,
    *,
    image_dimensions_xyz: Sequence[int],
    image_resolutions_mm: Sequence[float],
    template_resolutions_mm: Sequence[float],
    transform: ArrayLike,
    shape_slope: float = 2.0,
) -> NDArray[np.float32]:
    """Transform and normalize fixed rigid AtlasEM priors using the TOADS sampler.

    The selected ``DemonToadDeformableAtlas`` rigid branch uses the same
    precomputed transform, positive-support bounds, closest-edge trilinear
    sampler, and shape-slope restoration already ported for TOADS. This
    adapter computes those class bounds, reuses that source-checked sampler,
    then applies AtlasEM's masked prior normalization. It does not estimate
    the rigid transform or implement deformable registration.
    """
    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    active = np.asarray(mask, dtype=bool)
    if shapes.ndim != 4 or any(n < 2 for n in shapes.shape[:3]) or shapes.shape[3] < 1:
        raise ValueError("atlas_shapes must have shape (template_x, template_y, template_z, classes)")
    dimensions = tuple(int(value) for value in image_dimensions_xyz)
    if len(dimensions) != 3 or min(dimensions) < 2:
        raise ValueError("image_dimensions_xyz must contain three dimensions >= 2")
    if active.shape != dimensions:
        raise ValueError("mask must match image_dimensions_xyz")
    if not np.isfinite(shapes).all():
        raise ValueError("atlas_shapes must be finite")
    bounds = tuple(
        compute_toad_shape_bounds(
            shapes[..., k],
            image_dimensions_xyz=dimensions,
            image_resolutions_mm=tuple(image_resolutions_mm),
            template_resolutions_mm=tuple(template_resolutions_mm),
            transform=tuple(np.asarray(transform, dtype=np.float32)),
        )
        for k in range(shapes.shape[3])
    )
    transformed = transform_toads_shape_priors(
        shapes,
        image_dimensions_xyz=dimensions,
        image_resolutions_mm=tuple(image_resolutions_mm),
        template_resolutions_mm=tuple(template_resolutions_mm),
        transform=transform,
        class_bounds=bounds,
        shape_slope=shape_slope,
    )
    return normalize_atlas_em_priors(transformed, active)


def normalize_atlas_em_images(images: ArrayLike) -> AtlasEMNormalizedImages:
    """Apply ``AtlasEMPreprocess`` positive clipping and robust-span scaling.

    Source constructor maps negative values to zero, computes its robust min
    and max at ratio 0.01/refinement count 2, then divides by their difference
    without subtracting the low value. The input/output layout here is
    ``(x, y, z, channels)``; atlas mapping remains a separate step.
    """
    values = np.asarray(images, dtype=np.float32)
    if values.ndim != 4 or any(n < 1 for n in values.shape[:3]) or values.shape[3] < 1:
        raise ValueError("images must have shape (x, y, z, nonzero channels)")
    if not np.isfinite(values).all():
        raise ValueError("images must be finite")
    positive = np.maximum(values, np.float32(0.0)).astype(np.float32)
    channels = values.shape[3]
    low = np.empty(channels, np.float32)
    high = np.empty(channels, np.float32)
    normalized = np.empty_like(positive)
    for c in range(channels):
        low[c], high[c] = robust_intensity_bounds(positive[..., c], fraction=0.01, scales=2)
        span = np.float32(high[c] - low[c])
        if not np.isfinite(span) or span <= 0:
            raise ValueError(f"AtlasEM channel {c} has zero or invalid robust span")
        normalized[..., c] = np.asarray(positive[..., c] / span, dtype=np.float32)
    return AtlasEMNormalizedImages(normalized, low, high)


def initialize_atlas_em_centroids(
    robust_minima: ArrayLike,
    robust_maxima: ArrayLike,
    atlas_intensity_priors: ArrayLike,
) -> NDArray[np.float32]:
    """Apply active ``AtlasEMPreprocess.initialCentroids`` robustMaxMin formula."""
    low = np.asarray(robust_minima, dtype=np.float32)
    high = np.asarray(robust_maxima, dtype=np.float32)
    priors = np.asarray(atlas_intensity_priors, dtype=np.float32)
    if low.ndim != 1 or high.shape != low.shape or low.size < 1:
        raise ValueError("robust minima/maxima must be matching nonempty channel vectors")
    if priors.ndim != 2 or priors.shape[0] != low.size or priors.shape[1] < 1:
        raise ValueError("atlas intensity priors must have shape (channels, classes)")
    if not np.isfinite(low).all() or not np.isfinite(high).all() or not np.isfinite(priors).all():
        raise ValueError("robust bounds and atlas intensity priors must be finite")
    result = np.empty(priors.shape, dtype=np.float32)
    for c in range(low.size):
        span = np.float32(high[c] - low[c])
        if span <= 0:
            raise ValueError(f"AtlasEM channel {c} has a zero or negative robust span")
        offset = np.float32(low[c] / span)
        for k in range(priors.shape[1]):
            result[c, k] = np.float32(offset + priors[c, k])
    return result


def initialize_atlas_em_memberships(
    cropped_images: ArrayLike,
    centroids: ArrayLike,
    *,
    interior: Sequence[slice] | None = None,
) -> NDArray[np.float32]:
    """Reproduce source inverse-distance tissue memberships on a cropped grid.

    Voxels outside ``interior`` remain zero, matching the zero-initialized
    destination membership buffer outside the preprocessing crop. Exact
    centroid ties are grouped using the first channel's float equality, as in
    ``AtlasEMPreprocess.initialMemberships``.
    """
    values = np.asarray(cropped_images, dtype=np.float32)
    means = np.asarray(centroids, dtype=np.float32)
    if values.ndim != 4 or any(n < 2 for n in values.shape[:3]) or values.shape[3] < 1:
        raise ValueError("cropped_images must have shape (x, y, z, channels) with each spatial dimension >= 2")
    channels = values.shape[3]
    if means.ndim != 2 or means.shape[0] != channels or means.shape[1] < 1:
        raise ValueError("centroids must have shape (channels, classes)")
    if not np.isfinite(values).all() or not np.isfinite(means).all():
        raise ValueError("cropped images and centroids must be finite")
    # AtlasEMPreprocess uses xNi=nix-1 (and analogous y/z limits), exclusive.
    region = tuple(interior) if interior is not None else tuple(slice(0, n - 1) for n in values.shape[:3])
    if len(region) != 3 or any(s.start is None or s.stop is None or s.step not in (None, 1) or s.start < 0 or s.stop > size or s.start >= s.stop for s, size in zip(region, values.shape[:3], strict=True)):
        raise ValueError("interior must contain three nonempty unit-step slices within the image")
    classes = means.shape[1]
    group_sizes = np.zeros(classes, dtype=np.float32)
    for k in range(classes):
        for m in range(classes):
            if means[0, m] == means[0, k]:
                group_sizes[k] = np.float32(group_sizes[k] + np.float32(1.0))
    result = np.zeros((*values.shape[:3], classes), dtype=np.float32)
    zero = np.float32(1.0e-20)
    infinity = np.float32(1.0e20)
    for x in range(region[0].start, region[0].stop):
        for y in range(region[1].start, region[1].stop):
            for z in range(region[2].start, region[2].stop):
                denominator = np.float32(0.0)
                for k in range(classes):
                    distance = np.float32(0.0)
                    for c in range(channels):
                        difference = np.float32(values[x, y, z, c] - means[c, k])
                        distance = np.float32(distance + np.float32(difference * difference))
                    inverse_distance = np.float32(1.0) / distance if distance > zero else infinity
                    result[x, y, z, k] = inverse_distance
                    denominator = np.float32(denominator + inverse_distance)
                if denominator > 0:
                    for k in range(classes):
                        numerator = np.float32(group_sizes[k] * result[x, y, z, k])
                        result[x, y, z, k] = np.float32(numerator / denominator)
    return result


def initialize_atlas_em_memberships_from_descriptor(
    normalized: AtlasEMNormalizedImages,
    descriptor_path: str | Path,
    modalities: Sequence[str],
    *,
    interior: Sequence[slice] | None = None,
) -> AtlasEMDescriptorInitialization:
    """Build initial memberships from descriptor intensity priors in class order.

    The source `initialCentroids` combines per-modality robust bounds and the
    descriptor's intensity profile, then `initialMemberships` computes the
    inverse squared-distance memberships. The input is the output of
    :func:`normalize_atlas_em_images`; image loading, cropping, brain-mask
    handling, registration, and later EM updates remain separate stages.
    """
    from pymacruise.atlas_em_descriptor import parse_atlas_em_descriptor

    values = np.asarray(normalized.images, dtype=np.float32)
    names = tuple(str(modality) for modality in modalities)
    if values.ndim != 4 or len(names) != values.shape[3]:
        raise ValueError("modalities must name each normalized image channel in source order")
    aliases = {
        "T1": "T1_SPGR", "T1_SPGR": "T1_SPGR", "T2": "T2", "FLAIR": "FLAIR",
        "MPRAGE": "T1_MPRAGE", "T1_MPRAGE": "T1_MPRAGE",
        "MPRAGE_RAW": "T1_RAW", "T1_RAW": "T1_RAW", "PD": "PD",
        "PDFSE": "PDFSE", "DIR": "DIR",
    }
    descriptor = parse_atlas_em_descriptor(descriptor_path)
    try:
        profiles = np.asarray(
            [descriptor.intensity_atlas[aliases[name]] for name in names], dtype=np.float32
        )
    except KeyError as exc:
        raise ValueError(f"unsupported AtlasEM modality {exc.args[0]!r}") from exc
    centroids = initialize_atlas_em_centroids(
        normalized.robust_minima, normalized.robust_maxima, profiles
    )
    memberships = initialize_atlas_em_memberships(values, centroids, interior=interior)
    return AtlasEMDescriptorInitialization(tuple(item.name for item in descriptor.structures), centroids, memberships)
