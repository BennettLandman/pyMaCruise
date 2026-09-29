from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_data_terms import compute_lesion_toads_unfielded_data_energies


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_flair_weight/java_reference.json"


def _from_bits(values: list[int]) -> np.ndarray:
    return np.asarray(values, dtype=np.uint32).view(np.float32)


def _bits(values) -> list[int]:
    return [int(value) for value in np.asarray(values, dtype=np.float32).ravel().view(np.uint32)]


def test_flair_lesion_power_weight_matches_java_powertable() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert reference["powertable_sha256"] == "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
    assert len(reference["cases"]) == 6
    maximum_power_error = 0.0
    for case in reference["cases"]:
        result = compute_lesion_toads_unfielded_data_energies(
            _from_bits(case["image"]),
            _from_bits(case["class_centroid"]),
            _from_bits(case["lesion_centroid"]),
            np.zeros(3, dtype=np.float32),
            _from_bits(case["channel_weight"]),
            object_type="lesion",
            use_field=case["use_field"],
            intensity_field=_from_bits(case["field"]),
            flair_channel=int(case["flair_channel"]),
            use_lesion_weight=True,
            lesion_power_exponent=float(_from_bits([case["exponent_bits"]])[0]),
        )
        expected = _from_bits(case["energies"])
        assert _bits([result.class_energy]) == [int(case["energies"][0])]
        error = abs(float(result.lesion_energy) - float(expected[1]))
        maximum_power_error = max(maximum_power_error, error)
        assert error <= 1e-7, (case["sample"], case["use_field"], error)
        assert result.black_hole_energy == np.float32(0.0)
    assert maximum_power_error <= 1e-7
    assert any(case["squared_residual_bits"] > np.float32(4.0).view(np.uint32).item() for case in reference["cases"])


def test_flair_weight_requires_exponent_and_absent_flair_uses_unweighted_path() -> None:
    vector = np.asarray([0.7, 1.2], dtype=np.float32)
    centroid = np.asarray([0.3, 0.4], dtype=np.float32)
    lesion = np.asarray([0.2, 0.5], dtype=np.float32)
    weights = np.asarray([0.25, 0.75], dtype=np.float32)
    with pytest.raises(ValueError, match="lesion_power_exponent is required"):
        compute_lesion_toads_unfielded_data_energies(
            vector,
            centroid,
            lesion,
            np.zeros(2, dtype=np.float32),
            weights,
            object_type="lesion",
            flair_channel=1,
            use_lesion_weight=True,
        )
    fallback = compute_lesion_toads_unfielded_data_energies(
        vector,
        centroid,
        lesion,
        np.zeros(2, dtype=np.float32),
        weights,
        object_type="lesion",
        flair_channel=-1,
        use_lesion_weight=True,
    )
    ordinary = compute_lesion_toads_unfielded_data_energies(
        vector,
        centroid,
        lesion,
        np.zeros(2, dtype=np.float32),
        weights,
        object_type="lesion",
    )
    assert fallback == ordinary
