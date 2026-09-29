from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import convert_lesion_toads_class_membership


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_membership_conversion/java_reference.json"


def _from_bits(value: int) -> np.float32:
    return np.asarray(int(value) & 0xFFFFFFFF, dtype=np.uint32).view(np.float32)[()]


def _bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.int32))


def test_class_membership_conversion_matches_selected_java_float32() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert reference["numerics_sha256"] == "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb"
    assert len(reference["cases"]) == 11
    for case in reference["cases"]:
        result = convert_lesion_toads_class_membership(
            _from_bits(case["energy_bits"]),
            _from_bits(case["shape_bits"]),
            relation_type=case["relation"],
        )
        assert _bits(result) == case["membership_bits"], case

    cases = reference["cases"]
    assert cases[0]["membership_bits"] == cases[1]["membership_bits"]
    assert cases[2]["membership_bits"] != cases[1]["membership_bits"]
    assert cases[6]["membership_bits"] == cases[3]["membership_bits"]
    assert cases[7]["membership_bits"] == 0


@pytest.mark.parametrize(
    ("relation", "message"),
    [("faraway", "relation_type"), ("", "relation_type")],
)
def test_class_membership_conversion_rejects_unknown_relation(relation: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        convert_lesion_toads_class_membership(0.5, 0.2, relation_type=relation)
