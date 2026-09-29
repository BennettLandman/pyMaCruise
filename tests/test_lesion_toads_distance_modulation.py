from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import modulate_lesion_toads_membership_by_distances


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_distance_modulation/java_reference.json"


def _from_bits(value: int) -> np.float32:
    return np.asarray(value, dtype=np.uint32).view(np.float32)[()]


def _bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def test_lesion_distance_modulation_matches_selected_java_float32() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert len(reference["cases"]) == 12
    for case in reference["cases"]:
        result = modulate_lesion_toads_membership_by_distances(
            0.43,
            0.67,
            ventricle_distance=case["ventricle_distance"],
            gray_matter_distance=case["gray_matter_distance"],
            brainstem_distance=case["brainstem_distance"],
            max_ventricle_distance=case["max_ventricle_distance"],
            max_gray_matter_distance=case["max_gray_matter_distance"],
            max_brainstem_distance=case["max_brainstem_distance"],
            gm_flag=case["gm_flag"],
            spread=float(_from_bits(case["spread_bits"])),
        )
        assert _bits(result.membership) == case["membership_bits"], case
        assert _bits(result.relation) == case["relation_bits"], case


def test_lesion_distance_modulation_rejects_unsafe_distance_values() -> None:
    with pytest.raises(ValueError, match="nonnegative voxel distance"):
        modulate_lesion_toads_membership_by_distances(
            0.5, 1.0,
            ventricle_distance=-2,
            gray_matter_distance=-1,
            brainstem_distance=-1,
            max_ventricle_distance=9,
            max_gray_matter_distance=5,
            max_brainstem_distance=5,
            gm_flag=True,
            spread=0.5,
        )
