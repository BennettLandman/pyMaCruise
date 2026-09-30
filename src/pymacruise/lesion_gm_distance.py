"""Selected LesionTOADS gray-matter mask and capped byte distance support.

This ports the atlas-name selection in
``LesionToadSegmentation.computeGMDistance`` and reuses the selected byte
``fastDistance`` implementation. It is one support routine, not a full
LesionTOADS segmentation. Atlas labels are supplied in RAS-oriented Java
``(x, y, z)`` voxel order.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from pymacruise.lesion_distance import lesion_toads_distance_to_mask


_GM_DISTANCE_NAMES = frozenset(
    {
        "Cerebrum-GM",
        "CerebralGM",
        "CorticalGM",
        "Caudate",
        "GM",
        "Thalamus",
        "Putamen",
        "Sulcal-CSF",
        "SulcalCSF",
        "Subcortical-GM",
        "SubcorticalGM",
        "BasalGanglia",
        "GlobusPallidus",
    }
)


def lesion_toads_gm_distance(
    segmentation: np.ndarray,
    atlas_names: Sequence[str],
    max_distance: int,
) -> tuple[NDArray[np.bool_], NDArray[np.int8]]:
    """Build the source's exact GM-name mask and capped byte distance map.

    Matching is exact, case-sensitive, and punctuation-sensitive, as in
    Java ``String.equals``. ``segmentation`` contains integer indices into
    ``atlas_names`` and must already be RAS-oriented in Java ``(x, y, z)``
    order. Distance uses six-connected voxel steps; physical spacing does not
    enter the selected Java byte kernel. The return value retains its ``-1``
    unreachable sentinel and the Java kernel's `cap=1` extra-layer behavior.
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
        (index for index, name in enumerate(names) if name in _GM_DISTANCE_NAMES),
        dtype=np.intp,
    )
    mask = np.isin(labels, selected_ids)
    return mask, lesion_toads_distance_to_mask(mask, max_distance)
