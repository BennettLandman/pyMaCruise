from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_membership import update_lesion_toads_smoothing_energies


REFERENCE = (
    Path(__file__).parent
    / "data/synthetic/lesion_toads_smoothing_energy/java_reference.json"
)


def _from_bits(value: int) -> np.float32:
    return np.asarray(value & 0xFFFFFFFF, dtype=np.uint32).view(np.float32)[()]


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def test_six_face_smoothing_energy_matches_selected_java_float32() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert reference["face_order"] == ["X+", "X-", "Y+", "Y-", "Z+", "Z-"]
    assert len(reference["cases"]) == 4

    for case in reference["cases"]:
        result = update_lesion_toads_smoothing_energies(
            np.asarray([[_from_bits(v) for v in face] for face in case["class_neighbor_bits"]]),
            np.asarray([_from_bits(v) for v in case["lesion_neighbor_bits"]]),
            np.asarray([_from_bits(v) for v in case["black_hole_neighbor_bits"]]),
            class_index=case["class_index"],
            class_type=case["class_type"],
            smoothing_prior=float(_from_bits(case["smoothing_prior_bits"])),
            smoothing=float(_from_bits(case["smoothing_bits"])),
            class_energy=float(_from_bits(case["initial_class_energy_bits"])),
            lesion_energy=float(_from_bits(case["initial_lesion_energy_bits"])),
            black_hole_energy=float(_from_bits(case["initial_black_hole_energy_bits"])),
            outmask_membership_energy=float(_from_bits(case["initial_outmask_energy_bits"])),
        )
        assert _bits(result.class_energy) == case["class_energy_bits"], case
        assert _bits(result.lesion_energy) == case["lesion_energy_bits"], case
        assert _bits(result.black_hole_energy) == case["black_hole_energy_bits"], case
        assert _bits(result.outmask_membership_energy) == case["outmask_energy_bits"], case
        assert _bits(result.other_class_neighbor_squares) == case["other_class_neighbor_squares_bits"], case
        assert _bits(result.lesion_neighbor_squares) == case["lesion_neighbor_squares_bits"], case
        assert _bits(result.black_hole_neighbor_squares) == case["black_hole_neighbor_squares_bits"], case
        assert _bits(result.current_class_neighbor_squares) == case["current_class_neighbor_squares_bits"], case
