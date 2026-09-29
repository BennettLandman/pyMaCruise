from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import normalize_lesion_toads_voxel_memberships


REFERENCE = (
    Path(__file__).parent
    / "data/synthetic/lesion_toads_voxel_membership_normalization/java_reference.json"
)


def _from_bits(value: int) -> np.float32:
    return np.asarray(value & 0xFFFFFFFF, dtype=np.uint32).view(np.float32)[()]


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def test_voxel_membership_normalization_matches_selected_java_float32() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert len(reference["cases"]) == 8

    for case in reference["cases"]:
        result = normalize_lesion_toads_voxel_memberships(
            np.asarray([_from_bits(v) for v in case["class_membership_bits"]]),
            np.asarray([_from_bits(v) for v in case["class_energy_bits"]]),
            case["types"],
            lesion_membership=float(_from_bits(case["lesion_membership_bits"])),
            black_hole_membership=float(_from_bits(case["black_hole_membership_bits"])),
            outlier_membership=float(_from_bits(case["outlier_membership_bits"])),
            lesion_energy=float(_from_bits(case["lesion_energy_bits"])),
            black_hole_energy=float(_from_bits(case["black_hole_energy_bits"])),
        )
        assert _bits(result.denominator) == case["denominator_bits"], case
        assert [_bits(v) for v in result.class_memberships] == case["normalized_class_membership_bits"], case
        assert _bits(result.lesion_membership) == case["normalized_lesion_membership_bits"], case
        assert _bits(result.black_hole_membership) == case["normalized_black_hole_membership_bits"], case
        assert _bits(result.outlier_membership) == case["normalized_outlier_membership_bits"], case
        assert _bits(result.energy) == case["energy_bits"], case


@pytest.mark.parametrize(
    "memberships,energies,types",
    [
        (np.ones((1, 1)), np.ones(1), ["OBJECT"]),
        (np.ones(1), np.ones(2), ["OBJECT"]),
        (np.ones(1), np.ones(1), ["UNKNOWN"]),
        (np.ones(1), np.ones(1), []),
    ],
)
def test_voxel_membership_normalization_rejects_invalid_vectors(
    memberships: np.ndarray, energies: np.ndarray, types: list[str]
) -> None:
    with pytest.raises(ValueError):
        normalize_lesion_toads_voxel_memberships(
            memberships,
            energies,
            types,
            lesion_membership=0.0,
            black_hole_membership=0.0,
            outlier_membership=0.0,
            lesion_energy=0.0,
            black_hole_energy=0.0,
        )
