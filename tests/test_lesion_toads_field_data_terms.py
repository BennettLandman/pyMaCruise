from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_data_terms import compute_lesion_toads_unfielded_data_energies


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_field_data_terms/java_reference.json"


def _from_bits(values: list[int]) -> np.ndarray:
    return np.asarray(values, dtype=np.uint32).view(np.float32)


def _bits(values) -> list[int]:
    return [int(value) for value in np.asarray(values, dtype=np.float32).ravel().view(np.uint32)]


def test_field_corrected_lesion_data_energies_match_java_float32_formulas() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert len(reference["cases"]) == 9
    for case in reference["cases"]:
        result = compute_lesion_toads_unfielded_data_energies(
            _from_bits(case["image"]),
            _from_bits(case["class_centroid"]),
            _from_bits(case["lesion_centroid"]),
            _from_bits(case["black_hole_centroid"]),
            _from_bits(case["channel_weight"]),
            object_type=case["type"],
            use_field=True,
            intensity_field=_from_bits(case["field"]),
        )
        assert _bits(
            [result.class_energy, result.lesion_energy, result.black_hole_energy]
        ) == case["energies"], (case["type"], case["sample"])


def test_field_corrected_data_energies_require_matching_field() -> None:
    vector = np.ones(2, dtype=np.float32)
    with pytest.raises(ValueError, match="intensity_field is required"):
        compute_lesion_toads_unfielded_data_energies(
            vector,
            vector,
            vector,
            vector,
            vector,
            object_type="obj",
            use_field=True,
        )
