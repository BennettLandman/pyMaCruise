"""Source-compatible MIPAV binary morphology operations implemented with SciPy."""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import ndimage


def mipav_fill_holes_25d(image: ArrayLike) -> NDArray[np.bool_]:
    """Reproduce ``AlgorithmMorphology25D.FILL_HOLES`` for a binary 3D mask.

    MIPAV clears each 2D slice border to zero, flood-fills the 4-connected zero
    region from the upper-left seed, then maps reached voxels to zero and all
    remaining values to one. That border-clearing behavior differs from the
    usual ``scipy.ndimage.binary_fill_holes`` result, so it is preserved here.
    SciPy supplies only the per-slice 4-connected propagation primitive.
    """
    values = np.asarray(image)
    if values.ndim != 3 or any(size < 2 for size in values.shape):
        raise ValueError("image must be a 3D binary grid with each dimension >= 2")
    if not np.isfinite(values).all():
        raise ValueError("image must contain only finite values")
    foreground = values != 0
    background = ~foreground
    boundary = np.zeros(foreground.shape, dtype=bool)
    boundary[0, :, :] = True
    boundary[-1, :, :] = True
    boundary[:, 0, :] = True
    boundary[:, -1, :] = True
    # MIPAV explicitly sets every x/y slice boundary to the flood value 0.
    background[boundary] = True
    foreground[boundary] = False
    seed = np.zeros_like(background)
    seed[0, 0, :] = True
    cross = np.zeros((3, 3, 1), dtype=bool)
    cross[1, 1, 0] = True
    cross[0, 1, 0] = cross[2, 1, 0] = True
    cross[1, 0, 0] = cross[1, 2, 0] = True
    reachable = ndimage.binary_propagation(seed, structure=cross, mask=background)
    return np.asarray(~reachable, dtype=bool)
