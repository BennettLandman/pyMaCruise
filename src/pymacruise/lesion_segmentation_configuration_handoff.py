"""Active LesionTOADS segmentation configuration handoff.

Models the selected wrapper's immediate post-centroid calls that configure
lesion weighting, intensity ranges/scales, and gray-matter-distance usage.
The Java setters toggle one boolean or retain supplied float-array references;
they do not compute these inputs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass
class LesionToadsSegmentationConfiguration:
    """The four mutable fields affected by the wrapper handoff."""

    use_lesion_weight: bool = False
    intensity_max: NDArray[np.float32] | None = None
    intensity_scale: NDArray[np.float32] | None = None
    gm_flag: bool = False


def configure_lesion_toads_segmentation(
    state: LesionToadsSegmentationConfiguration,
    *,
    lesion_weight: bool,
    intensity_max: NDArray[np.float32],
    intensity_scale: NDArray[np.float32],
) -> None:
    """Apply the wrapper's four source setters in their call order.

    ``intensity_max`` and ``intensity_scale`` are retained by reference, just
    like Java ``setIntensityMax`` and ``setIntensityScale``. ``gm_flag`` is
    always enabled on this wrapper path. The existing ``lesion_state``
    constructor model contains initial ones for the intensity arrays and a
    channel-weight vector, but not these final references or either boolean.
    """
    if lesion_weight:
        state.use_lesion_weight = True
    else:
        state.use_lesion_weight = False
    state.intensity_max = intensity_max
    state.intensity_scale = intensity_scale
    state.gm_flag = True
