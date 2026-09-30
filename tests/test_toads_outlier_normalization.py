"""Selected-TOADS source check for reuse of the existing voxel normalizer."""

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_membership import normalize_lesion_toads_voxel_memberships


REFERENCE = Path(__file__).parent / "data/synthetic/toads_outlier_normalization/java_reference.json"


def _values(bits: list[int]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32)


def _bits(value: float | np.ndarray) -> list[int]:
    return np.asarray(value, dtype=np.float32).ravel().view(np.int32).tolist()


def test_selected_toads_outlier_normalization_matches_existing_shared_helper() -> None:
    reference = json.loads(REFERENCE.read_text())
    for case in reference["cases"]:
        result = normalize_lesion_toads_voxel_memberships(
            _values(case["membership_bits"]),
            _values(case["energy_bits"]),
            case["class_types"],
            lesion_membership=0.0,
            black_hole_membership=0.0,
            outlier_membership=float(_values([case["outlier_bits"]])[0]),
            lesion_energy=0.0,
            black_hole_energy=0.0,
        )
        assert _bits(result.denominator) == [case["denominator_bits"]]
        assert _bits(result.class_memberships) == case["normalized_bits"]
        assert _bits(result.outlier_membership) == [case["normalized_outlier_bits"]]
        assert _bits(result.energy) == [case["objective_bits"]]
