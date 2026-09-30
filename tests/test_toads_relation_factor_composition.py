"""Selected Java check for composed TOADS per-class relation weights."""

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_relations import toads_membership_relation_factors


REFERENCE = Path(__file__).parent / "data/synthetic/toads_relation_factor_composition/java_reference.json"


def _float_values(bits: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32).reshape(shape)


def _bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).ravel().view(np.int32).tolist()


def test_composed_relation_weights_match_selected_java_bits() -> None:
    reference = json.loads(REFERENCE.read_text())
    shape = tuple(reference["shape"])
    classes = reference["classes"]
    closest = np.asarray(reference["closest_labels"], dtype=np.int8).reshape((*shape, 4))
    transformed = _float_values(reference["shape_bits"], (*shape, classes))
    closest_copy = closest.copy()
    transformed_copy = transformed.copy()

    for case in reference["priors"]:
        prior = _float_values([case["prior_bits"]], ())[()]
        result = toads_membership_relation_factors(closest, transformed, float(prior))
        assert _bits(result) == case["factor_bits"]
        assert result.dtype == np.float32
    np.testing.assert_array_equal(closest, closest_copy)
    np.testing.assert_array_equal(transformed, transformed_copy)


def test_composed_relation_weights_reject_invalid_prior() -> None:
    closest = np.full((2, 2, 2, 4), -1, dtype=np.int8)
    shapes = np.zeros((2, 2, 2, 3), dtype=np.float32)
    with pytest.raises(ValueError, match="prior_factor"):
        toads_membership_relation_factors(closest, shapes, float("nan"))
