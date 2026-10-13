from __future__ import annotations

import runpy
from pathlib import Path


FIXTURE_GENERATOR = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "plus_atlas_consistent"
    / "generate_fixture.py"
)


def test_atlas_topology_to_maseg_mapping_preserves_cortical_tissue_classes() -> None:
    mapping = runpy.run_path(str(FIXTURE_GENERATOR))["TOPOLOGY_ID_TO_MASEG"]

    assert mapping == {
        1: 0,
        5: 5,
        6: 5,
        14: 75,
        15: 100,
        16: 75,
        17: 75,
        18: 75,
        22: 75,
        24: 75,
        25: 44,
        26: 44,
        27: 44,
    }
    assert [key for key, value in mapping.items() if value == 100] == [15]
    assert {key for key, value in mapping.items() if value == 44} == {25, 26, 27}
