"""Selected LesionTOADS ventricle mask and capped distance support.

This ports the exact atlas-name selection in the historically misspelled
``LesionToadSegmentation.computeVentirclsDistance`` method and delegates its
byte distance output to the selected ``fastDistance`` implementation. It is
not a LesionTOADS segmentation implementation. Arrays use RAS-oriented Java
``(x, y, z)`` voxel order.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from pymacruise.lesion_distance import lesion_toads_distance_to_mask


_VENTRICLE_NAMES = frozenset(
    {"Ventricles", "LeftVentricle", "RightVentricle", "ThirdVentricle"}
)


def lesion_toads_ventricle_distance(
    segmentation: np.ndarray,
    atlas_names: Sequence[str],
    max_distance: int,
) -> tuple[NDArray[np.bool_], NDArray[np.int8]]:
    """Build the selected ventricle mask and capped byte distance map.

    A voxel is selected only when its atlas name exactly matches one of
    ``Ventricles``, ``LeftVentricle``, ``RightVentricle``, or
    ``ThirdVentricle``. Matching is case- and punctuation-sensitive, matching
    Java ``String.equals``. ``segmentation`` contains integer indices into
    ``atlas_names`` and must already be in the source's RAS-oriented
    ``(x, y, z)`` voxel order.

    Distances use six-connected voxel steps, ignore physical spacing, and
    retain the source byte kernel's ``-1`` sentinel and cap behavior,
    including its extra propagation when ``max_distance == 1``.
    """
    labels = np.asarray(segmentation)
    if labels.ndim != 3 or any(size == 0 for size in labels.shape):
        raise ValueError("segmentation must be a nonempty 3D (x, y, z) array")
    if not np.issubdtype(labels.dtype, np.integer):
        raise TypeError("segmentation must contain integer atlas indices")

    names = tuple(atlas_names)
    if any(not isinstance(name, str) for name in names):
        raise TypeError("atlas_names must contain only strings")
    if labels.size and (int(labels.min()) < 0 or int(labels.max()) >= len(names)):
        raise ValueError("segmentation contains an index outside atlas_names")

    selected_ids = np.fromiter(
        (index for index, name in enumerate(names) if name in _VENTRICLE_NAMES),
        dtype=np.intp,
    )
    mask = np.isin(labels, selected_ids)
    return mask, lesion_toads_distance_to_mask(mask, max_distance)
