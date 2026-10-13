"""Compiled exact-order registered atlas-prior sampling for LesionTOADS."""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True, fastmath=False, nogil=True)
def registered_shapes_numba(
    atlas_maps: np.ndarray,
    matrix: np.ndarray,
    bounds: np.ndarray,
    image_dimensions_xyz: tuple[int, int, int],
) -> np.ndarray:
    """Apply the selected scalar transform/interpolation in Java loop order.

    Float32 rounding, inclusive interpolation bounds, per-class support bounds,
    interpolation term order, and class-0/outside defaults match the Python
    scalar oracle and the selected ImageFunctions operation.
    """
    nx, ny, nz = image_dimensions_xyz
    class_count = atlas_maps.shape[0]
    output = np.empty((nx, ny, nz, class_count), dtype=np.float32)
    for x in range(nx):
        x_float = np.float32(x)
        for y in range(ny):
            y_float = np.float32(y)
            for z in range(nz):
                z_float = np.float32(z)
                mapped_x = np.float32(matrix[0, 0] * x_float)
                mapped_x = np.float32(mapped_x + np.float32(matrix[0, 1] * y_float))
                mapped_x = np.float32(mapped_x + np.float32(matrix[0, 2] * z_float))
                mapped_x = np.float32(mapped_x + matrix[0, 3])
                mapped_y = np.float32(matrix[1, 0] * x_float)
                mapped_y = np.float32(mapped_y + np.float32(matrix[1, 1] * y_float))
                mapped_y = np.float32(mapped_y + np.float32(matrix[1, 2] * z_float))
                mapped_y = np.float32(mapped_y + matrix[1, 3])
                mapped_z = np.float32(matrix[2, 0] * x_float)
                mapped_z = np.float32(mapped_z + np.float32(matrix[2, 1] * y_float))
                mapped_z = np.float32(mapped_z + np.float32(matrix[2, 2] * z_float))
                mapped_z = np.float32(mapped_z + matrix[2, 3])

                for k in range(class_count):
                    default = np.float32(1.0 if k == 0 else 0.0)
                    if (
                        x < bounds[k, 0] or x >= bounds[k, 3]
                        or y < bounds[k, 1] or y >= bounds[k, 4]
                        or z < bounds[k, 2] or z >= bounds[k, 5]
                    ):
                        output[x, y, z, k] = default
                        continue
                    if (
                        mapped_x < np.float32(0.0)
                        or mapped_x > np.float32(atlas_maps.shape[1] - 2)
                        or mapped_y < np.float32(0.0)
                        or mapped_y > np.float32(atlas_maps.shape[2] - 2)
                        or mapped_z < np.float32(0.0)
                        or mapped_z > np.float32(atlas_maps.shape[3] - 2)
                    ):
                        output[x, y, z, k] = default
                        continue

                    ix = int(np.floor(mapped_x))
                    iy = int(np.floor(mapped_y))
                    iz = int(np.floor(mapped_z))
                    a = np.float32(mapped_x - np.float32(ix))
                    b = np.float32(mapped_y - np.float32(iy))
                    c = np.float32(mapped_z - np.float32(iz))
                    na = np.float32(np.float32(1.0) - a)
                    nb = np.float32(np.float32(1.0) - b)
                    nc = np.float32(np.float32(1.0) - c)
                    image = atlas_maps[k]
                    terms = np.empty(8, dtype=np.float32)
                    terms[0] = np.float32(np.float32(np.float32(na * nb) * nc) * image[ix, iy, iz])
                    terms[1] = np.float32(np.float32(np.float32(a * nb) * nc) * image[ix + 1, iy, iz])
                    terms[2] = np.float32(np.float32(np.float32(na * b) * nc) * image[ix, iy + 1, iz])
                    terms[3] = np.float32(np.float32(np.float32(na * nb) * c) * image[ix, iy, iz + 1])
                    terms[4] = np.float32(np.float32(np.float32(a * b) * nc) * image[ix + 1, iy + 1, iz])
                    terms[5] = np.float32(np.float32(np.float32(na * b) * c) * image[ix, iy + 1, iz + 1])
                    terms[6] = np.float32(np.float32(np.float32(a * nb) * c) * image[ix + 1, iy, iz + 1])
                    terms[7] = np.float32(np.float32(np.float32(a * b) * c) * image[ix + 1, iy + 1, iz + 1])
                    result = np.float32(terms[0] + terms[1])
                    for term in range(2, 8):
                        result = np.float32(result + terms[term])
                    output[x, y, z, k] = result
    return output
