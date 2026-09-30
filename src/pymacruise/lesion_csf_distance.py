"""Selected LesionTOADS CSF-name mask and capped distance support.

This ports the dormant ``LesionToadSegmentation.computeCSFDistance`` helper
from the selected VandyCRUISE source. The selected wrapper comments out its
call, initializes ``CSFDist`` to zero, and does not export this distance map;
this module therefore adds no active workflow behavior. Distance uses the
existing byte ``fastDistance`` port and RAS-oriented Java ``(x, y, z)`` voxel
order.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from pymacruise.lesion_distance import lesion_toads_distance_to_mask


_CSF_DISTANCE_NAMES = frozenset(
    {
        "Caudate",
        "Thalamus",
        "Putamen",
        "Subcortical-GM",
        "SubcorticalGM",
        "BasalGanglia",
    }
)


def lesion_toads_csf_distance(
    segmentation: np.ndarray,
    atlas_names: Sequence[str],
    max_distance: int = 5,
) -> tuple[NDArray[np.bool_], NDArray[np.int8]]:
    """Return the exact selected CSF seed mask and its byte distance field.

    ``segmentation`` contains integer indices into ``atlas_names`` and must
    already be in RAS-oriented Java ``(x, y, z)`` order. Names use exact,
    case-sensitive Java ``String.equals`` matching. If the name table has
    duplicate selected names, every corresponding index is a seed, as the
    source checks the name at each voxel. The selected field starts with a
    cap of 5; callers may provide the configured cap explicitly. The returned
    distance retains the shared kernel's signed-byte ``-1`` sentinel and
    cap-one extra propagation behavior. Physical voxel spacing is ignored.
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
        (index for index, name in enumerate(names) if name in _CSF_DISTANCE_NAMES),
        dtype=np.intp,
    )
    mask = np.isin(labels, selected_ids)
    return mask, lesion_toads_distance_to_mask(mask, max_distance)
