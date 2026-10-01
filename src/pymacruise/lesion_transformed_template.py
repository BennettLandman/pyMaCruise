"""Rigid LesionTOADS atlas-label resampling onto the current image grid.

Ports the selected ``LesionToadAtlas.computeTransformedTemplate`` method for
single rigid transforms. Labels are sampled with the source's bounded nearest
neighbor rule and ``label[0]`` background fill. The active
``AlgorithmLesionToads`` wrapper calls it after its optional rigid
``alignShapes`` registration, recrop, and ``setImageInfo`` sequence.

Arrays are RAS ``[x,y,z]`` and resolutions are positive mm/voxel. This helper
does not handle multi-transform averaging or non-rigid transform modes.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_crop import toads_rotation_matrix


def transform_lesion_toads_template(
    template: ArrayLike,
    transform: ArrayLike,
    image_shape: tuple[int, int, int],
    *,
    image_resolution_mm: tuple[float, float, float],
    template_resolution_mm: tuple[float, float, float],
    background_label: int,
) -> NDArray[Any]:
    """Sample atlas labels on the post-registration image grid.

    `transform` is the six float32 rigid parameter vector used by VandyCRUISE:
    inverse-Rodrigues parameters followed by physical translations. Centering,
    axis order, float32 arithmetic, and half-voxel rounding follow
    ``ParametricTransform.imageToTemplate`` and
    ``ImageFunctions.nearestNeighborInterpolation``. ``background_label`` is
    the atlas metadata's first class value (Java ``label[0]``), not the label
    at atlas voxel ``[0,0,0]``. The transform mode is intentionally restricted
    to rigid, which is what the selected wrapper initializes and updates for
    LesionTOADS shape registration.
    """
    atlas = np.asarray(template)
    trans = np.asarray(transform, dtype=np.float32)
    if atlas.ndim != 3 or any(n < 1 for n in atlas.shape):
        raise ValueError("template must be a nonempty 3-D atlas label array")
    if trans.shape != (6,) or not np.isfinite(trans).all():
        raise ValueError("transform must contain six finite float32 values")
    if len(image_shape) != 3 or any(
        isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) or n < 1
        for n in image_shape
    ):
        raise ValueError("image_shape must contain three positive integer dimensions")
    image_shape = tuple(int(n) for n in image_shape)
    image_res = np.asarray(image_resolution_mm, dtype=np.float32)
    template_res = np.asarray(template_resolution_mm, dtype=np.float32)
    for name, values in (("image_resolution_mm", image_res), ("template_resolution_mm", template_res)):
        if values.shape != (3,) or not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError(f"{name} must contain three finite positive values")
    if not isinstance(background_label, (int, np.integer)) or not -128 <= int(background_label) <= 127:
        raise ValueError("background_label must fit a Java signed byte")
    if atlas.dtype.kind not in "iu":
        raise ValueError("template must contain integer label values")

    rotation = toads_rotation_matrix(trans[:3])
    image_center = tuple(np.float32(size / 2.0) for size in image_shape)
    template_center = tuple(np.float32(size / 2.0) for size in atlas.shape)
    # The Java atlas stores its topology as byte[][][]; return int8 to preserve
    # signed-byte labels while keeping all supported values exact.
    output = np.full(image_shape, np.int8(background_label), dtype=np.int8)
    atlas_i8 = atlas.astype(np.int8, copy=False)
    nx, ny, nz = image_shape
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                coords: list[np.float32] = []
                # ParametricTransform.imageToTemplate evaluates Xi/Yi/Zi
                # before rotation; preserve Java's float32 operation order.
                source_delta = (
                    np.float32(np.float32(np.float32(x) - image_center[0]) * image_res[0]),
                    np.float32(np.float32(np.float32(y) - image_center[1]) * image_res[1]),
                    np.float32(np.float32(np.float32(z) - image_center[2]) * image_res[2]),
                )
                for axis in range(3):
                    total = np.float32(0.0)
                    for component in range(3):
                        product = np.float32(rotation[axis, component] * source_delta[component])
                        total = np.float32(total + product)
                    total = np.float32(total + trans[axis + 3])
                    total = np.float32(total / template_res[axis])
                    coords.append(np.float32(total + template_center[axis]))

                # ImageFunctions rejects continuous out-of-bounds coordinates
                # before rounding, then Numerics.round uses ties away from 0.
                if any(coords[a] < 0 or coords[a] > atlas.shape[a] - 1 for a in range(3)):
                    continue
                indices = tuple(
                    int(np.floor(c + np.float32(0.5))) if c > 0
                    else int(np.ceil(c - np.float32(0.5)))
                    for c in coords
                )
                output[x, y, z] = atlas_i8[indices]
    return output
