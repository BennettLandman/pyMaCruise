"""Source-convention crop and topology sampling for the selected TOADS path.

These helpers operate on arrays already normalized to the caller's RAS grid.
The ``orientation_signs`` argument preserves the source MIPAV axis-direction
factors; it is not a NIfTI orientation detector. Coordinates and dimensions use
the selected Java ``[x][y][z]`` convention, with channels leading for images.

Paper context: crop/pad and topology-template resampling are setup operations
for the atlas-driven multiple-object workflow discussed by Bazin & Pham (2008).
Their concrete bounds and coordinate conventions follow the selected Java
wrapper and atlas source, not a standalone prescription in the paper. See
``docs/papers/toads.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ToadsCropResult:
    """Padded image, inclusive source bounds, and post-crop rigid transform."""

    images: np.ndarray
    bounds: tuple[int, int, int, int, int, int]
    transform: np.ndarray


def toads_rotation_matrix(parameters: np.ndarray) -> np.ndarray:
    """Recreate ``RotationMatrix.setParameters`` (inverse Rodrigues form)."""
    values = np.asarray(parameters, dtype=np.float32)
    if values.shape != (3,) or not np.isfinite(values).all():
        raise ValueError("TOADS rotation parameters must be three finite values")
    a, b, c = (float(value) for value in values)
    norm = float(np.sqrt(a * a + b * b + c * c))
    if norm > 1.0:
        magnitude = norm % 2.0
        scale = (magnitude - 2.0 if magnitude > 1.0 else magnitude) / norm
        a, b, c = a * scale, b * scale, c * scale
    a2, b2, c2 = a * a, b * b, c * c
    w_sq = 1.0 - a2 - b2 - c2
    if w_sq < 0 and w_sq > -1e-10:
        w_sq = 0.0
    if w_sq < 0:
        raise ValueError("TOADS rotation parameters produce an invalid rotation")
    w = float(np.sqrt(w_sq))
    return np.asarray(
        [
            [1.0 - 2.0 * b2 - 2.0 * c2, 2.0 * a * b - 2.0 * w * c, 2.0 * a * c + 2.0 * w * b],
            [2.0 * a * b + 2.0 * w * c, 1.0 - 2.0 * a2 - 2.0 * c2, 2.0 * b * c - 2.0 * w * a],
            [2.0 * a * c - 2.0 * w * b, 2.0 * b * c + 2.0 * w * a, 1.0 - 2.0 * a2 - 2.0 * b2],
        ],
        dtype=np.float32,
    )


def find_toads_crop_bounds(
    images: np.ndarray,
    topology: np.ndarray,
    transform: np.ndarray,
    *,
    image_resolution: tuple[float, float, float],
    template_resolution: tuple[float, float, float],
    orientation_signs: tuple[float, float, float],
    robust_ranges: np.ndarray,
    background_ratio: float,
    background_label: int = 1,
    normalized: bool = True,
) -> tuple[int, int, int, int, int, int]:
    """Find source-compatible inclusive image bounds from image and atlas data."""
    image = np.asarray(images, dtype=np.float32)
    atlas = np.asarray(topology)
    trans = np.asarray(transform, dtype=np.float32)
    ranges = np.asarray(robust_ranges, dtype=np.float32)
    image_res = np.asarray(image_resolution, dtype=np.float32)
    template_res = np.asarray(template_resolution, dtype=np.float32)
    signs = np.asarray(orientation_signs, dtype=np.float32)
    if image.ndim != 4 or any(size < 1 for size in image.shape[1:]):
        raise ValueError("TOADS crop images must have shape (channels, x, y, z)")
    if atlas.ndim != 3:
        raise ValueError("TOADS topology must be a 3D [x,y,z] array")
    if trans.shape != (6,) or ranges.shape != (image.shape[0], 2):
        raise ValueError("TOADS transform or robust ranges have invalid shape")
    for name, values in (("images", image), ("transform", trans), ("ranges", ranges)):
        if not np.isfinite(values).all():
            raise ValueError(f"TOADS {name} contain NaN or infinite values")
    for name, values in (("image_resolution", image_res), ("template_resolution", template_res)):
        if values.shape != (3,) or not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError(f"TOADS {name} must be three finite positive values")
    if signs.shape != (3,) or not np.isfinite(signs).all() or np.any(np.abs(signs) != 1):
        raise ValueError("TOADS orientation_signs must contain three values from {-1, 1}")
    if not np.isfinite(background_ratio) or background_ratio < 0:
        raise ValueError("background_ratio must be finite and nonnegative")

    nx, ny, nz = image.shape[1:]
    bounds = [nx - 1, 0, ny - 1, 0, nz - 1, 0]
    robust_spans = ranges[:, 1] - ranges[:, 0]
    threshold = background_ratio * robust_spans
    for channel in range(image.shape[0]):
        include = (image[channel] > background_ratio) if normalized else False
        include = np.logical_or(include, image[channel] > threshold[channel])
        points = np.argwhere(include)
        if points.size:
            for axis in range(3):
                bounds[2 * axis] = min(bounds[2 * axis], int(points[:, axis].min()))
                bounds[2 * axis + 1] = max(bounds[2 * axis + 1], int(points[:, axis].max()))

    rot = toads_rotation_matrix(trans[:3])
    x_center, y_center, z_center = (np.float32(size / 2.0) for size in atlas.shape)
    image_centers = (np.float32(nx / 2.0), np.float32(ny / 2.0), np.float32(nz / 2.0))
    for x, y, z in np.argwhere(atlas > background_label):
        atlas_delta = np.asarray(
            [
                np.float32((x - x_center) * template_res[0] - trans[3]),
                np.float32((y - y_center) * template_res[1] - trans[4]),
                np.float32((z - z_center) * template_res[2] - trans[5]),
            ], dtype=np.float32,
        )
        # Java's computeImageCoordinates multiplies by rotation columns here.
        mapped = np.empty(3, dtype=np.float32)
        for axis in range(3):
            total = np.float32(0.0)
            for row in range(3):
                total = np.float32(total + np.float32(rot[row, axis] * atlas_delta[row]))
            mapped[axis] = np.float32(total * signs[axis] / image_res[axis] + image_centers[axis])
        for axis, coord in enumerate(mapped):
            low_index, high_index = 2 * axis, 2 * axis + 1
            low = int(np.floor(coord))
            high = int(np.ceil(coord))
            bounds[low_index] = min(bounds[low_index], low)
            bounds[high_index] = max(bounds[high_index], high)

    if any(bounds[2 * axis] > bounds[2 * axis + 1] for axis in range(3)):
        raise ValueError("TOADS crop has no image or topology support")
    return tuple(bounds)


def crop_toads_images(
    images: np.ndarray,
    bounds: tuple[int, int, int, int, int, int],
    padding: int,
    transform: np.ndarray,
    *,
    image_resolution: tuple[float, float, float],
    orientation_signs: tuple[float, float, float],
) -> ToadsCropResult:
    """Crop and zero-pad channels, then update the source rigid translation."""
    source = np.asarray(images, dtype=np.float32)
    if source.ndim != 4:
        raise ValueError("TOADS crop images must have shape (channels, x, y, z)")
    if len(bounds) != 6 or padding < 0:
        raise ValueError("TOADS crop bounds or padding are invalid")
    x0, x_n, y0, y_n, z0, z_n = (int(v) for v in bounds)
    if x0 > x_n or y0 > y_n or z0 > z_n:
        raise ValueError("TOADS crop bounds must be inclusive and ordered")
    trans = np.array(transform, dtype=np.float32, copy=True)
    if trans.shape != (6,) or not np.isfinite(trans).all():
        raise ValueError("TOADS transform must contain six finite values")
    image_res = np.asarray(image_resolution, dtype=np.float32)
    signs = np.asarray(orientation_signs, dtype=np.float32)
    if image_res.shape != (3,) or np.any(image_res <= 0) or signs.shape != (3,):
        raise ValueError("TOADS crop geometry resolutions or orientation are invalid")
    nx, ny, nz = source.shape[1:]
    shape = (x_n - x0 + 1 + 2 * padding, y_n - y0 + 1 + 2 * padding, z_n - z0 + 1 + 2 * padding)
    rot = toads_rotation_matrix(trans[:3])
    axis_delta = np.asarray(
        [
            (shape[0] / 2.0 + x0 - padding - nx / 2.0) * image_res[0] / signs[0],
            (shape[1] / 2.0 + y0 - padding - ny / 2.0) * image_res[1] / signs[1],
            (shape[2] / 2.0 + z0 - padding - nz / 2.0) * image_res[2] / signs[2],
        ], dtype=np.float32,
    )
    for axis in range(3):
        delta = np.float32(0.0)
        for component in range(3):
            delta = np.float32(delta + np.float32(rot[axis, component] * axis_delta[component]))
        trans[axis + 3] = np.float32(trans[axis + 3] + delta)

    cropped = np.zeros((source.shape[0], *shape), dtype=np.float32)
    for x in range(x0, x_n + 1):
        for y in range(y0, y_n + 1):
            for z in range(z0, z_n + 1):
                if 0 <= x < nx and 0 <= y < ny and 0 <= z < nz:
                    cropped[:, x - x0 + padding, y - y0 + padding, z - z0 + padding] = source[:, x, y, z]
    return ToadsCropResult(cropped, tuple(bounds), trans)


def uncrop_toads_images(
    cropped: np.ndarray,
    original_shape: tuple[int, int, int],
    bounds: tuple[int, int, int, int, int, int],
    padding: int,
    transform: np.ndarray,
    *,
    image_resolution: tuple[float, float, float],
    orientation_signs: tuple[float, float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """Reinsert the padded crop into a zero image and undo its translation."""
    source = np.asarray(cropped, dtype=np.float32)
    if source.ndim != 4 or len(original_shape) != 3:
        raise ValueError("TOADS uncrop expects channel-first and 3D shapes")
    x0, x_n, y0, y_n, z0, z_n = bounds
    nx, ny, nz = (int(v) for v in original_shape)
    expected = (x_n - x0 + 1 + 2 * padding, y_n - y0 + 1 + 2 * padding, z_n - z0 + 1 + 2 * padding)
    if source.shape[1:] != expected:
        raise ValueError("TOADS uncrop crop shape does not match its bounds and padding")
    trans = np.array(transform, dtype=np.float32, copy=True)
    if trans.shape != (6,) or not np.isfinite(trans).all():
        raise ValueError("TOADS transform must contain six finite values")
    image_res = np.asarray(image_resolution, dtype=np.float32)
    signs = np.asarray(orientation_signs, dtype=np.float32)
    shape = source.shape[1:]
    rot = toads_rotation_matrix(trans[:3])
    axis_delta = np.asarray(
        [
            (shape[0] / 2.0 + x0 - padding - nx / 2.0) * image_res[0] / signs[0],
            (shape[1] / 2.0 + y0 - padding - ny / 2.0) * image_res[1] / signs[1],
            (shape[2] / 2.0 + z0 - padding - nz / 2.0) * image_res[2] / signs[2],
        ], dtype=np.float32,
    )
    for axis in range(3):
        delta = np.float32(0.0)
        for component in range(3):
            delta = np.float32(delta + np.float32(rot[axis, component] * axis_delta[component]))
        trans[axis + 3] = np.float32(trans[axis + 3] - delta)

    output = np.zeros((source.shape[0], nx, ny, nz), dtype=np.float32)
    starts = (max(x0 - padding, 0), max(y0 - padding, 0), max(z0 - padding, 0))
    ends = (min(x_n + padding, nx - 1), min(y_n + padding, ny - 1), min(z_n + padding, nz - 1))
    if all(start <= end for start, end in zip(starts, ends)):
        for x in range(starts[0], ends[0] + 1):
            for y in range(starts[1], ends[1] + 1):
                for z in range(starts[2], ends[2] + 1):
                    output[:, x, y, z] = source[:, x - x0 + padding, y - y0 + padding, z - z0 + padding]
    return output, trans


def uncrop_toads_buffer(
    cropped: np.ndarray,
    original_shape: tuple[int, int, int],
    bounds: tuple[int, int, int, int, int, int],
    padding: int,
    background: float | int,
) -> np.ndarray:
    """Return a selected-Java ``uncropAndBuffer`` buffer.

    A 3D ``[x,y,z]`` input yields one scalar volume; a 4D ``[channel,x,y,z]``
    input yields channel-major volumes. In each volume x is the fastest axis,
    matching ``x + nx*y + nx*ny*z`` in ``MultipleToadPreprocess``. The returned
    array is one-dimensional and uses the input dtype, so byte labels retain
    Java byte bit patterns and float values retain their source precision.
    """
    source = np.asarray(cropped)
    if source.ndim not in (3, 4) or not np.issubdtype(source.dtype, np.number):
        raise ValueError("TOADS buffer input must be a numeric [x,y,z] or [channel,x,y,z] array")
    if len(original_shape) != 3 or len(bounds) != 6:
        raise ValueError("TOADS buffer shape and bounds must describe three spatial axes")
    if padding < 0:
        raise ValueError("TOADS buffer padding must be nonnegative")
    nx, ny, nz = (int(v) for v in original_shape)
    x0, x_n, y0, y_n, z0, z_n = (int(v) for v in bounds)
    if min(nx, ny, nz) < 1 or x0 > x_n or y0 > y_n or z0 > z_n:
        raise ValueError("TOADS buffer source geometry is invalid")
    expected = (x_n - x0 + 1 + 2 * padding, y_n - y0 + 1 + 2 * padding, z_n - z0 + 1 + 2 * padding)
    spatial = source.shape[-3:]
    if spatial != expected:
        raise ValueError("TOADS buffer crop shape does not match its bounds and padding")

    channels = 1 if source.ndim == 3 else source.shape[0]
    if channels < 1:
        raise ValueError("TOADS buffer must contain at least one channel")
    output = np.full((channels, nx, ny, nz), background, dtype=source.dtype)
    starts = (max(x0 - padding, 0), max(y0 - padding, 0), max(z0 - padding, 0))
    ends = (min(x_n + padding, nx - 1), min(y_n + padding, ny - 1), min(z_n + padding, nz - 1))
    if all(start <= end for start, end in zip(starts, ends)):
        cropped_starts = (starts[0] - x0 + padding, starts[1] - y0 + padding, starts[2] - z0 + padding)
        cropped_ends = (ends[0] - x0 + padding + 1, ends[1] - y0 + padding + 1, ends[2] - z0 + padding + 1)
        slices = tuple(slice(start, end) for start, end in zip(starts, tuple(end + 1 for end in ends)))
        crop_slices = tuple(slice(start, end) for start, end in zip(cropped_starts, cropped_ends))
        if source.ndim == 3:
            output[0][slices] = source[crop_slices]
        else:
            output[(slice(None), *slices)] = source[(slice(None), *crop_slices)]
    # C-order flattening of [channel,z,y,x] produces Java's channel-major,
    # x-fastest linear layout while preserving the [x,y,z] public convention.
    if source.ndim == 3:
        return output[0].transpose(2, 1, 0).reshape(-1)
    return output.transpose(0, 3, 2, 1).reshape(-1)


def transform_toads_topology(
    topology: np.ndarray,
    bounds: tuple[int, int, int, int, int, int],
    padding: int,
    transform: np.ndarray,
    *,
    crop_shape: tuple[int, int, int],
    image_resolution: tuple[float, float, float],
    template_resolution: tuple[float, float, float],
    orientation_signs: tuple[float, float, float],
    background_label: int = 1,
) -> np.ndarray:
    """Nearest-neighbor resample categorical atlas labels onto padded crop."""
    atlas = np.asarray(topology)
    if atlas.ndim != 3:
        raise ValueError("TOADS topology must be 3D")
    x0, x_n, y0, y_n, z0, z_n = bounds
    expected = (x_n - x0 + 1 + 2 * padding, y_n - y0 + 1 + 2 * padding, z_n - z0 + 1 + 2 * padding)
    if tuple(crop_shape) != expected:
        raise ValueError("TOADS transformed topology crop shape does not match bounds")
    trans = np.asarray(transform, dtype=np.float32)
    if trans.shape != (6,) or not np.isfinite(trans).all():
        raise ValueError("TOADS transform must contain six finite values")
    image_res = np.asarray(image_resolution, dtype=np.float32)
    template_res = np.asarray(template_resolution, dtype=np.float32)
    signs = np.asarray(orientation_signs, dtype=np.float32)
    rot = toads_rotation_matrix(trans[:3])
    centers_image = tuple(np.float32(dimension / 2.0) for dimension in crop_shape)
    centers_template = tuple(np.float32(dimension / 2.0) for dimension in atlas.shape)
    output = np.full(expected, background_label, dtype=atlas.dtype)
    for x in range(x0, x_n + 1):
        for y in range(y0, y_n + 1):
            for z in range(z0, z_n + 1):
                source_coords = (x, y, z)
                target = []
                for axis in range(3):
                    total = np.float32(0.0)
                    for component in range(3):
                        term = np.float32(rot[axis, component] * np.float32(source_coords[component] - centers_image[component]))
                        term = np.float32(term * image_res[component] / signs[component])
                        total = np.float32(total + term)
                    total = np.float32(total + trans[axis + 3])
                    target.append(np.float32(total / template_res[axis] + centers_template[axis]))
                if any(target[axis] < 0 or target[axis] > atlas.shape[axis] - 1 for axis in range(3)):
                    value = background_label
                else:
                    indices = [int(np.floor(v + np.float32(0.5))) if v > 0 else int(np.ceil(v - np.float32(0.5))) for v in target]
                    value = atlas[indices[0], indices[1], indices[2]]
                output[x - x0 + padding, y - y0 + padding, z - z0 + padding] = value
    return output
