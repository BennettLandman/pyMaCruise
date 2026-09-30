"""Selected-Java comparison for a composed TOADS OUTLIER voxel update."""

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_math import toads_outlier_membership_update


REFERENCE = Path(__file__).parent / "data/synthetic/toads_outlier_membership_update/java_reference.json"


def test_outlier_voxel_update_matches_selected_java_bits() -> None:
    reference = json.loads(REFERENCE.read_text())
    for case in reference["cases"]:
        v = np.asarray(case["input_bits"], dtype=np.int32).view(np.float32)
        class_value, outlier_value, energy = toads_outlier_membership_update(
            v[0], v[1], v[2], v[3],
            smoothing_prior=v[4], smoothing=v[5], neighbor_sum=v[6],
            prior_coefficient=v[7], inverse_prior=v[8],
        )
        assert int(class_value.view(np.int32)) == case["class_bits"]
        assert int(outlier_value.view(np.int32)) == case["outlier_bits"]
        assert int(energy.view(np.int32)) == case["energy_bits"]


def test_outlier_voxel_update_rejects_nonfinite_values() -> None:
    with pytest.raises(ValueError, match="inputs must be finite"):
        toads_outlier_membership_update(0.5, 0.7, 0.25, 1.0, neighbor_sum=np.nan)
