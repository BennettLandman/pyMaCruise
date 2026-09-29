from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_membership import update_lesion_toads_outmask_membership


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_outmask_membership/java_reference.json"


def _from_bits(value: int) -> np.float32:
    return np.asarray(value & 0xFFFFFFFF, dtype=np.uint32).view(np.float32)[()]


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def test_outmask_membership_matches_selected_java_float32() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert len(reference["cases"]) == 12
    for case in reference["cases"]:
        result = update_lesion_toads_outmask_membership(
            float(_from_bits(case["class_input_bits"])),
            float(_from_bits(case["outmask_input_bits"])),
            float(_from_bits(case["outlier_input_bits"])),
            relation_type=case["relation_type"],
            shape_prior=float(_from_bits(case["shape_bits"])),
        )
        assert _bits(result.class_membership) == case["class_output_bits"], case
        assert _bits(result.outlier_membership) == case["outlier_output_bits"], case
