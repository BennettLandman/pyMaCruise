"""Selected-Java comparison for the TOADS OUTLIER membership mix."""

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_math import toads_outlier_membership_mix


REFERENCE = Path(__file__).parent / "data/synthetic/toads_outlier_membership_mix/java_reference.json"


def test_outlier_membership_mix_matches_selected_java_bits() -> None:
    reference = json.loads(REFERENCE.read_text())
    for case in reference["cases"]:
        inputs = np.asarray(case["input_bits"], dtype=np.int32).view(np.float32)
        class_value, outlier_value = toads_outlier_membership_mix(*inputs)
        assert int(class_value.view(np.int32)) == case["class_bits"]
        assert int(outlier_value.view(np.int32)) == case["outlier_bits"]


def test_outlier_membership_mix_rejects_nonfinite_inputs() -> None:
    with pytest.raises(ValueError, match="inputs must be finite"):
        toads_outlier_membership_mix(0.5, 0.25, np.inf)
