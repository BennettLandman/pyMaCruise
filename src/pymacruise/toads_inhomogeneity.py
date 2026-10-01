"""Source-ordered Chebyshev basis used by TOADS' polynomial bias field.

The selected `InhomogeneityCorrection2` builds this 3D polynomial at each
voxel before assembling a global weighted least-squares field. This module
ports only that basis primitive; it does not fit or apply a correction field.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def toads_chebyshev_basis(
    voxel_xyz: ArrayLike,
    image_dimensions_xyz: tuple[int, int, int],
    *,
    degree: int,
) -> NDArray[np.float32]:
    """Return the selected Java polynomial basis for one integer voxel.

    Degree 1 has four terms, degree 2 has ten, and degree 3 has twenty. The
    source normalizes each coordinate as `index / (dimension - 1)` (so the
    domain is [0, 1], not the canonical [-1, 1]) and stores mixed terms in its
    explicit historical order. The wrapper's selected maximum degree is 3.
    """
    point = np.asarray(voxel_xyz)
    if point.shape != (3,) or not np.issubdtype(point.dtype, np.integer):
        raise ValueError("voxel_xyz must contain three integer indices")
    dims = np.asarray(image_dimensions_xyz)
    if dims.shape != (3,) or not np.issubdtype(dims.dtype, np.integer) or np.any(dims < 2):
        raise ValueError("image_dimensions_xyz must contain three integer sizes >= 2")
    if np.any(point < 0) or np.any(point >= dims):
        raise ValueError("voxel_xyz must lie within image_dimensions_xyz")
    if isinstance(degree, bool) or degree not in (1, 2, 3):
        raise ValueError("degree must be one of the selected supported values: 1, 2, or 3")

    x, y, z = (np.float32(point[a]) / np.float32(dims[a] - 1) for a in range(3))
    terms = np.zeros(20, dtype=np.float32)
    terms[0] = np.float32(1.0)
    terms[1:4] = (x, y, z)
    if degree >= 2:
        terms[4] = np.float32(np.float32(2.0) * x * x - np.float32(1.0))
        terms[5] = np.float32(np.float32(2.0) * y * y - np.float32(1.0))
        terms[6] = np.float32(np.float32(2.0) * z * z - np.float32(1.0))
        terms[7] = np.float32(terms[1] * terms[2])
        terms[8] = np.float32(terms[2] * terms[3])
        terms[9] = np.float32(terms[3] * terms[1])
    if degree >= 3:
        terms[10] = np.float32(np.float32(np.float32(4.0) * x * x * x) - np.float32(3.0) * x)
        terms[11] = np.float32(np.float32(np.float32(4.0) * y * y * y) - np.float32(3.0) * y)
        terms[12] = np.float32(np.float32(np.float32(4.0) * z * z * z) - np.float32(3.0) * z)
        terms[13] = np.float32(terms[4] * terms[2])
        terms[14] = np.float32(terms[1] * terms[5])
        terms[15] = np.float32(terms[5] * terms[3])
        terms[16] = np.float32(terms[2] * terms[6])
        terms[17] = np.float32(terms[6] * terms[1])
        terms[18] = np.float32(terms[3] * terms[4])
        terms[19] = np.float32(np.float32(terms[1] * terms[2]) * terms[3])
    return terms[: (4, 10, 20)[degree - 1]].copy()
