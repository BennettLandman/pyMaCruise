from __future__ import annotations

import numpy as np

from pymacruise.tgdm_diagnostic import _array_comparison, levelset_component_summary


def test_tgdm_diagnostic_counts_connectivity_dependent_islands() -> None:
    level_set = np.ones((9, 9, 9), dtype=np.float32)
    level_set[2, 2, 2] = -1.0
    level_set[3, 3, 3] = -1.0
    level_set[5:7, 5:7, 5:7] = -1.0

    summary = levelset_component_summary(level_set)

    assert summary["negative_voxels"] == 10
    assert summary["connectivity_6"]["components"] == 3
    assert summary["connectivity_6"]["singleton_components"] == 2
    assert summary["connectivity_18"]["components"] == 3
    assert summary["connectivity_26"]["components"] == 2
    assert summary["connectivity_26"]["singleton_components"] == 0
    assert summary["connectivity_26"]["largest_component_voxels"] == 8


def test_tgdm_diagnostic_separates_legacy_uncrop_high_faces() -> None:
    replay = np.zeros((5, 6, 7), dtype=np.float32)
    saved = replay.copy()
    saved[-1, :, :] = 20
    saved[:, -1, :] = 20
    saved[:, :, -1] = 20

    comparison = _array_comparison(replay, saved)

    assert comparison["legacy_uncrop_written_region"]["bitwise_equal"]
    assert comparison["legacy_uncrop_written_region"]["different_voxels"] == 0
    assert comparison["exclusive_high_boundary_planes"]["different_voxels"] > 0
