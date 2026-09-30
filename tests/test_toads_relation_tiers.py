"""Selected Java comparison for TOADS class support-tier assignment."""

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_relations import toads_membership_relation_tiers


REFERENCE = Path(__file__).parent / "data/synthetic/toads_relation_tiers/java_reference.json"


def _float_values(bits: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32).reshape(shape)


def test_membership_relation_tiers_match_selected_java() -> None:
    reference = json.loads(REFERENCE.read_text())
    shape = tuple(reference["shape"])
    classes = reference["classes"]
    closest = np.asarray(reference["closest_labels"], dtype=np.int8).reshape((*shape, 4))
    transformed = _float_values(reference["shape_bits"], (*shape, classes))
    closest_copy = closest.copy()
    transformed_copy = transformed.copy()

    relation_types, eligible = toads_membership_relation_tiers(closest, transformed)

    np.testing.assert_array_equal(relation_types.ravel(), reference["relation_types"])
    np.testing.assert_array_equal(eligible.ravel().astype(np.int8), reference["eligible"])
    assert relation_types.dtype == np.int8
    assert eligible.dtype == np.bool_
    np.testing.assert_array_equal(closest, closest_copy)
    np.testing.assert_array_equal(transformed, transformed_copy)


def test_membership_relation_tiers_reject_malformed_inputs() -> None:
    closest = np.full((2, 2, 2, 4), -1, dtype=np.int8)
    shape = np.ones((2, 2, 2, 3), dtype=np.float32)
    with pytest.raises(ValueError, match="four slots"):
        toads_membership_relation_tiers(closest[..., :3], shape)
    with pytest.raises(ValueError, match="match the closest-label grid"):
        toads_membership_relation_tiers(closest, np.ones((3, 2, 2, 3), dtype=np.float32))
    with pytest.raises(ValueError, match="valid class indices"):
        toads_membership_relation_tiers(np.full_like(closest, 3), shape)
    invalid_shapes = shape.copy()
    invalid_shapes[0, 0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        toads_membership_relation_tiers(closest, invalid_shapes)
