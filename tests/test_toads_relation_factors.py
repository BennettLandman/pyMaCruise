"""Selected Java float32 comparisons for TOADS relation factors."""

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_relations import (
    toads_distance_factor,
    toads_neighbor_factor,
    toads_shape_factor,
)


REFERENCE = Path(__file__).parent / "data/synthetic/toads_relation_factors/java_reference.json"


def _float_from_bits(bits: int) -> np.float32:
    return np.asarray([bits], dtype=np.int32).view(np.float32)[0]


def _bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).ravel().view(np.int32).tolist()


def test_distance_factors_match_selected_java_float_bits() -> None:
    reference = json.loads(REFERENCE.read_text())
    for case in reference["distances"]:
        limit = _float_from_bits(case["limit_bits"])
        result = toads_distance_factor(np.asarray(case["input"], dtype=np.int8), float(limit))
        assert _bits(result) == case["factor_bits"]


def test_shape_and_neighbor_factors_match_selected_java_float_bits() -> None:
    reference = json.loads(REFERENCE.read_text())
    for case in reference["shapes"]:
        shape = _float_from_bits(case["shape_bits"])
        prior = _float_from_bits(case["prior_bits"])
        assert _bits(toads_shape_factor(shape, float(prior))) == [case["shape_factor_bits"]]
        assert _bits(toads_neighbor_factor(shape, float(prior))) == [case["neighbor_factor_bits"]]


@pytest.mark.parametrize("limit", [0.0, -1.0, float("nan"), float("inf")])
def test_distance_factor_rejects_nonpositive_or_nonfinite_limit(limit: float) -> None:
    with pytest.raises(ValueError, match="marching_limit"):
        toads_distance_factor(np.asarray([0, 1], dtype=np.int8), limit)


def test_factors_reject_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="signed-byte"):
        toads_distance_factor(np.asarray([128], dtype=np.int16), 1.0)
    with pytest.raises(ValueError, match="integer"):
        toads_distance_factor(np.asarray([0.0]), 1.0)
    with pytest.raises(ValueError, match="finite"):
        toads_shape_factor(np.asarray([np.nan], dtype=np.float32), 1.0)
    with pytest.raises(ValueError, match="finite"):
        toads_neighbor_factor(np.asarray([0.5], dtype=np.float32), float("inf"))
