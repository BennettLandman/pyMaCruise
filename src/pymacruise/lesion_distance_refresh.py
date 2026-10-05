"""Refresh LesionTOADS spatial distance state from its current segmentation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from .lesion_brainstem import lesion_toads_brainstem_distance
from .lesion_gm_distance import lesion_toads_gm_distance
from .lesion_registered_initial_membership import LesionInitialDistanceMaps
from .lesion_ventricle import lesion_toads_ventricle_distance


@dataclass(frozen=True)
class LesionDistanceRefresh:
    """The maps rebuilt by one GM-enabled source membership pass."""

    maps: LesionInitialDistanceMaps
    stage_order: tuple[str, ...]


def refresh_lesion_toads_distance_maps(
    segmentation: NDArray,
    atlas_names: Sequence[str],
    *,
    maximum_ventricle_distance: int,
    maximum_gm_distance: int,
    maximum_brainstem_distance: int,
    previous_maps: LesionInitialDistanceMaps | None = None,
) -> LesionDistanceRefresh:
    """Rebuild exact selected VandyCRUISE masks/maps in source call order.

    Java ``computeMemberships`` refreshes ventricle, brainstem, then GM maps
    from the current segmentation before its GM-enabled membership update.
    Each adapter preserves the selected exact-name and byte-distance semantics.
    """
    labels = np.asarray(segmentation)
    if previous_maps is None and min(
        maximum_ventricle_distance, maximum_gm_distance, maximum_brainstem_distance
    ) == 0:
        raise ValueError("previous_maps are required when a source distance cap disables refresh")
    refresh_order: list[str] = []
    if maximum_ventricle_distance != 0:
        ventricle_mask, ventricle = lesion_toads_ventricle_distance(
            labels, atlas_names, maximum_ventricle_distance
        )
        refresh_order.append("computeVentirclsDistance")
    else:
        assert previous_maps is not None
        ventricle_mask, ventricle = previous_maps.ventricle_mask, previous_maps.ventricle
    if maximum_brainstem_distance != 0:
        brainstem_mask, brainstem = lesion_toads_brainstem_distance(
            labels, atlas_names, maximum_brainstem_distance
        )
        refresh_order.append("computeBstemDistance")
    else:
        assert previous_maps is not None
        brainstem_mask, brainstem = previous_maps.brainstem_mask, previous_maps.brainstem
    if maximum_gm_distance != 0:
        gm_mask, gray_matter = lesion_toads_gm_distance(
            labels, atlas_names, maximum_gm_distance
        )
        refresh_order.append("computeGMDistance")
    else:
        assert previous_maps is not None
        gm_mask, gray_matter = previous_maps.gray_matter_mask, previous_maps.gray_matter
    return LesionDistanceRefresh(
        LesionInitialDistanceMaps(
            ventricle_mask, ventricle, gm_mask, gray_matter,
            brainstem_mask, brainstem,
        ),
        tuple(refresh_order),
    )
