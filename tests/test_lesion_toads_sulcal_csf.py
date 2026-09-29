from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_data_terms import compute_lesion_toads_unfielded_data_energies


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_sulcal_csf/java_reference.json"


def _from_bits(values: list[int]) -> np.ndarray:
    return np.asarray(values, dtype=np.uint32).view(np.float32)


def test_sulcal_csf_mprage_energy_matches_selected_java() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert reference["powertable_sha256"] == "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
    assert len(reference["cases"]) == 9
    for case in reference["cases"]:
        result = compute_lesion_toads_unfielded_data_energies(
            _from_bits(case["image"]),
            _from_bits(case["class_centroid"]),
            np.zeros(3, dtype=np.float32),
            np.zeros(3, dtype=np.float32),
            _from_bits(case["channel_weight"]),
            object_type="obj",
            use_field=case["use_field"],
            intensity_field=_from_bits(case["field"]),
            structure_name=case["name"],
            t1_channel=int(case["t1_channel"]),
            has_mprage=case["has_mprage"],
            csf_power_exponent=float(_from_bits([case["exponent_bits"]])[0]),
        )
        assert int(np.asarray(result.class_energy, dtype=np.float32).view(np.uint32)) == case["energy_bits"], case
        assert result.lesion_energy == np.float32(0.0)
        assert result.black_hole_energy == np.float32(0.0)


def test_sulcal_csf_mprage_requires_power_exponent_and_valid_t1_index() -> None:
    vector = np.ones(2, dtype=np.float32)
    with pytest.raises(ValueError, match="csf_power_exponent is required"):
        compute_lesion_toads_unfielded_data_energies(
            vector,
            vector,
            vector,
            vector,
            vector,
            object_type="obj",
            structure_name="Sulcal-CSF",
            t1_channel=0,
            has_mprage=True,
        )
    with pytest.raises(ValueError, match="t1_channel"):
        compute_lesion_toads_unfielded_data_energies(
            vector,
            vector,
            vector,
            vector,
            vector,
            object_type="obj",
            t1_channel=2,
        )
