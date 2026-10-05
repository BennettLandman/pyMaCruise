from __future__ import annotations

import numpy as np

from pymacruise.lesion_brainstem import lesion_toads_brainstem_distance
from pymacruise.lesion_distance_refresh import refresh_lesion_toads_distance_maps
from pymacruise.lesion_gm_distance import lesion_toads_gm_distance
from pymacruise.lesion_registered_initial_membership import LesionInitialDistanceMaps
from pymacruise.lesion_ventricle import lesion_toads_ventricle_distance


def test_gm_distance_refresh_rebuilds_source_masks_in_call_order() -> None:
    names = ("Background", "Ventricles", "Cerebrum-GM", "Brainstem")
    segmentation = np.zeros((5, 5, 5), dtype=np.int16)
    segmentation[1, 2, 2] = 1
    segmentation[2, 2, 2] = 2
    segmentation[3, 2, 2] = 3
    previous = LesionInitialDistanceMaps(
        np.zeros_like(segmentation, dtype=bool),
        np.full_like(segmentation, -1, dtype=np.int8),
        np.zeros_like(segmentation, dtype=bool),
        np.full_like(segmentation, -1, dtype=np.int8),
        np.zeros_like(segmentation, dtype=bool),
        np.full_like(segmentation, -1, dtype=np.int8),
    )

    refreshed = refresh_lesion_toads_distance_maps(
        segmentation,
        names,
        maximum_ventricle_distance=2,
        maximum_gm_distance=3,
        maximum_brainstem_distance=25,
        previous_maps=previous,
    )

    assert refreshed.stage_order == (
        "computeVentirclsDistance", "computeBstemDistance", "computeGMDistance"
    )
    for mask, distance, adapter, cap in (
        (refreshed.maps.ventricle_mask, refreshed.maps.ventricle, lesion_toads_ventricle_distance, 2),
        (refreshed.maps.brainstem_mask, refreshed.maps.brainstem, lesion_toads_brainstem_distance, 25),
        (refreshed.maps.gray_matter_mask, refreshed.maps.gray_matter, lesion_toads_gm_distance, 3),
    ):
        expected_mask, expected_distance = adapter(segmentation, names, cap)
        np.testing.assert_array_equal(mask, expected_mask)
        np.testing.assert_array_equal(distance, expected_distance)
        assert distance.dtype == np.int8


def test_zero_source_distance_caps_keep_the_previous_map_state() -> None:
    segmentation = np.zeros((3, 3, 3), dtype=np.int16)
    previous = LesionInitialDistanceMaps(
        np.ones_like(segmentation, dtype=bool),
        np.full_like(segmentation, 1, dtype=np.int8),
        np.ones_like(segmentation, dtype=bool),
        np.full_like(segmentation, 2, dtype=np.int8),
        np.ones_like(segmentation, dtype=bool),
        np.full_like(segmentation, 3, dtype=np.int8),
    )

    refreshed = refresh_lesion_toads_distance_maps(
        segmentation,
        ("Background",),
        maximum_ventricle_distance=0,
        maximum_gm_distance=0,
        maximum_brainstem_distance=0,
        previous_maps=previous,
    )

    assert refreshed.stage_order == ()
    for name in (
        "ventricle_mask", "ventricle", "gray_matter_mask", "gray_matter",
        "brainstem_mask", "brainstem",
    ):
        np.testing.assert_array_equal(getattr(refreshed.maps, name), getattr(previous, name))
