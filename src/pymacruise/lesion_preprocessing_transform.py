"""Selected LesionTOADS fixed-transform image resampling helper.

This ports the dormant ``LesionToadPreprocess.transformNewImage`` utility.
It samples a template-grid image onto the inclusive, padded image-space crop
using the supplied rigid transform. The selected wrapper has no call site for
this helper, so this module exposes its numerical operation without claiming
that it runs in the active LesionTOADS workflow.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.lesion_shapes import _linear_interpolation_with_default
from pymacruise.toads_crop import toads_rotation_matrix


def transform_lesion_toads_image(
    template_image: ArrayLike,
    transform: ArrayLike,
    image_shape: tuple[int, int, int],
    crop_bounds: tuple[int, int, int, int, int, int],
    *,
    image_resolution: tuple[float, float, float],
    template_resolution: tuple[float, float, float],
    orientation_signs: tuple[float, float, float] = (1.0, 1.0, 1.0),
    padding: int = 0,
    background: float = 0.0,
) -> NDArray[np.float32]:
    """Resample the template-grid image using source float32/order semantics.

    ``template_image`` is a 3-D ``(x, y, z)`` RAS array. ``image_shape`` and
    inclusive ``crop_bounds=(x0,xN,y0,yN,z0,zN)`` describe the target image
    grid. The six transform parameters are inverse-Rodrigues rotations then
    millimetre translations. Image/template centers are derived as each
    dimension divided by two, matching ``LesionToadPreprocess`` construction.
    Output shape includes ``padding`` voxels on each side; those voxels retain
    ``background``. Sampling outside the valid linear-interpolation domain
    also returns ``background``.

    This helper takes arrays and explicit geometry; callers should load/store
    NIfTI through ``image_io`` so spatial images remain compressed and RAS.
    """
    source = np.asarray(template_image, dtype=np.float32)
    trans = np.asarray(transform, dtype=np.float32)
    if source.ndim != 3 or any(n < 2 for n in source.shape):
        raise ValueError("template_image must be a 3-D array of size at least two per axis")
    if trans.shape != (6,) or not np.isfinite(trans).all():
        raise ValueError("transform must contain six finite float32 values")
    if len(image_shape) != 3 or any(not isinstance(n, (int, np.integer)) or n < 1 for n in image_shape):
        raise ValueError("image_shape must contain three positive integer dimensions")
    if len(crop_bounds) != 6:
        raise ValueError("crop_bounds must contain inclusive low/high indices for three axes")
    bounds = tuple(int(n) for n in crop_bounds)
    shape = tuple(int(n) for n in image_shape)
    for axis in range(3):
        low, high = bounds[2 * axis : 2 * axis + 2]
        if low < 0 or high < low or high >= shape[axis]:
            raise ValueError("crop_bounds must be ordered and inside image_shape")
    if isinstance(padding, (bool, np.bool_)) or not isinstance(padding, (int, np.integer)) or padding < 0:
        raise ValueError("padding must be a nonnegative integer")
    if not np.isfinite(background):
        raise ValueError("background must be finite")

    image_res = np.asarray(image_resolution, dtype=np.float32)
    template_res = np.asarray(template_resolution, dtype=np.float32)
    signs = np.asarray(orientation_signs, dtype=np.float32)
    for name, values in (("image_resolution", image_res), ("template_resolution", template_res)):
        if values.shape != (3,) or not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError(f"{name} must contain three finite positive values")
    if signs.shape != (3,) or not np.isfinite(signs).all() or np.any(np.abs(signs) != 1):
        raise ValueError("orientation_signs must contain three values from {-1, 1}")

    rotation = toads_rotation_matrix(trans[:3])
    image_center = tuple(np.float32(n / 2.0) for n in shape)
    template_center = tuple(np.float32(n / 2.0) for n in source.shape)
    lows = bounds[::2]
    highs = bounds[1::2]
    pad = int(padding)
    output_shape = tuple(highs[a] - lows[a] + 1 + 2 * pad for a in range(3))
    output = np.full(output_shape, np.float32(background), dtype=np.float32)

    for x in range(lows[0], highs[0] + 1):
        for y in range(lows[1], highs[1] + 1):
            for z in range(lows[2], highs[2] + 1):
                target = (x, y, z)
                template_coordinates: list[np.float32] = []
                for axis in range(3):
                    total = np.float32(0.0)
                    for component in range(3):
                        delta = np.float32(
                            np.float32(target[component]) - image_center[component]
                        )
                        product = np.float32(rotation[axis, component] * delta)
                        product = np.float32(product * image_res[component])
                        product = np.float32(product / signs[component])
                        total = np.float32(total + product)
                    total = np.float32(total + trans[axis + 3])
                    total = np.float32(total / template_res[axis])
                    template_coordinates.append(
                        np.float32(total + template_center[axis])
                    )

                sampled = _linear_interpolation_with_default(
                    source,
                    tuple(template_coordinates),
                    np.float32(background),
                )
                output[x - lows[0] + pad, y - lows[1] + pad, z - lows[2] + pad] = sampled
    return output
