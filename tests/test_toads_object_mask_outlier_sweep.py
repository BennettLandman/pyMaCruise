"""Selected Java comparison for the composed TOADS OBJECT/OUTLIER sweep."""

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_math import toads_object_mask_outlier_membership_sweep


REFERENCE = Path(__file__).parent / "data/synthetic/toads_object_mask_outlier_sweep/java_reference.json"


def _values(bits: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32).reshape(shape)


def _bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).ravel().view(np.int32).tolist()


def test_object_outlier_membership_sweep_matches_selected_java() -> None:
    reference = json.loads(REFERENCE.read_text())
    spatial = tuple(reference["shape"])
    classes = reference["classes"]
    class_shape = (*spatial, classes)
    result, outlier, objective = toads_object_mask_outlier_membership_sweep(
        _values(reference["class_data_energy_bits"], class_shape),
        _values(reference["outlier_data_energy_bits"], spatial),
        _values(reference["initial_memberships_bits"], class_shape),
        _values(reference["initial_outlier_membership_bits"], spatial),
        reference["class_types"],
        _values(reference["shape_priors_bits"], class_shape),
        np.asarray(reference["closest_labels"], dtype=np.int8).reshape((*spatial, 4)),
        _values(reference["smoothing_prior_bits"], (classes,)),
        reference["smoothing"],
        reference["prior_factor"],
    )
    assert _bits(result) == reference["updated_bits"]
    assert _bits(outlier) == reference["updated_outlier_bits"]
    assert _bits(np.asarray([objective])) == [reference["objective_bits"]]
    assert result.dtype == outlier.dtype == np.float32

