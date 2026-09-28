"""Source-faithful distance fields used by the MaCRUISE dura stage.

The selected `MedicAlgorithmRemoveDura09` implementation delegates signed
distance construction to JIST `DistanceField`. The implementation is shared
with the verified six-neighbor fast-marching port used elsewhere in the
package; Java stores every distance update as float32, so this wrapper selects
that arithmetic explicitly. Distances and caps are in voxel units, matching
the source class (which does not use image voxel spacing).
"""

from __future__ import annotations

import numpy as np

from pymacruise.ace import _signed_distance_fmm


def dura_signed_distance_field(
    level_set: np.ndarray,
    max_distance_voxels: float,
) -> np.ndarray:
    """Return the selected JIST signed-distance field for a 3D level set.

    Inputs are expected in the already RAS-normalized array order used by the
    workflow. The algorithm operates on voxel indices and deliberately does
    not apply physical image spacing. Values are cast to float32, and the
    result is a float32 array with the source class's finite-distance cap.
    """
    return _signed_distance_fmm(
        level_set,
        max_distance=max_distance_voxels,
        round_updates_to_float32=True,
    )
