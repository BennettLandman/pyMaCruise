"""Selected LesionTOADS brainstem-related mask and capped distance support.

This ports only ``LesionToadSegmentation.computeBstemDistance``'s exact atlas
name selection and delegates the distance map to the selected byte-valued
``fastDistance`` helper. It is not a LesionTOADS segmentation implementation.
Arrays are in RAS voxel order and ``(x, y, z)`` index order, as in the source.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from pymacruise.lesion_distance import lesion_toads_distance_to_mask


_BRAINSTEM_DISTANCE_NAMES = frozenset(
    {"Brainstem", "InterVentricular-WM", "FornixWM"}
)


def lesion_toads_brainstem_distance(
    segmentation: np.ndarray,
    atlas_names: Sequence[str],
    max_distance: int,
) -> tuple[NDArray[np.bool_], NDArray[np.int8]]:
    """Build the selected brainstem-related mask and capped byte distance map.

    A voxel belongs to the mask only when its atlas name is exactly
    ``Brainstem``, ``InterVentricular-WM``, or ``FornixWM``. Names are not
    normalized and aliases are not accepted, matching Java ``String.equals``.
    ``segmentation`` contains integer indices into ``atlas_names`` and must
    already be RAS-oriented in Java-style ``(x, y, z)`` order.

    The returned distance map uses six-connected voxel steps, byte dtype, and
    the source's ``-1`` sentinel and cap behavior. Physical spacing is ignored
    by the selected Java kernel.
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
        (index for index, name in enumerate(names) if name in _BRAINSTEM_DISTANCE_NAMES),
        dtype=np.intp,
    )
    mask = np.isin(labels, selected_ids)
    return mask, lesion_toads_distance_to_mask(mask, max_distance)
