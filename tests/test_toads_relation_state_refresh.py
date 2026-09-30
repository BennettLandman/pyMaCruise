"""Source-ordered synthetic checks for composed TOADS relation state."""

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_relations import toads_refresh_relation_state


REFERENCE = Path(__file__).parent / "data/synthetic/toads_relation_state_refresh/java_reference.json"


def test_relation_state_refresh_matches_selected_java_composition() -> None:
    reference = json.loads(REFERENCE.read_text())
    shape = tuple(reference["shape"])
    classes = reference["classes"]
    segmentation = np.asarray(reference["segmentation"], dtype=np.int8).reshape(shape)
    membership_bits = np.asarray(reference["membership_float_bits"], dtype=np.uint32)
    memberships = membership_bits.view(np.float32).reshape((*shape, classes))
    initial = np.asarray(reference["initial_boundary_distances"], dtype=np.int8).reshape((*shape, classes))
    initial_copy = initial.copy()

    result = toads_refresh_relation_state(
        segmentation, memberships, initial, reference["limit"]
    )

    np.testing.assert_array_equal(result.relation.relation_map.ravel(), reference["relation_map"])
    np.testing.assert_array_equal(
        result.relation.relation_membership.ravel().astype(np.int8),
        reference["relation_membership"],
    )
    np.testing.assert_array_equal(result.boundary_distances.ravel(), reference["boundary_distances"])
    np.testing.assert_array_equal(result.closest_labels.ravel(), reference["closest_labels"])
    assert result.relation.diffusion_iterations == reference["relation_diffusion_iterations"]
    assert result.boundary_diffusion_passes == reference["boundary_diffusion_passes"] == 1
    assert result.closest_label_passes == 5 * reference["limit"] - 1
    assert result.relation.relation_map.dtype == np.int16
    assert result.boundary_distances.dtype == np.int8
    assert result.closest_labels.dtype == np.int8
    np.testing.assert_array_equal(initial, initial_copy)


@pytest.mark.parametrize("limit", [-1, 1.5, True])
def test_relation_state_refresh_rejects_invalid_limit(limit: object) -> None:
    segmentation = np.zeros((3, 3, 3), dtype=np.int8)
    memberships = np.ones((3, 3, 3, 2), dtype=np.float32)
    boundary = np.zeros_like(memberships, dtype=np.int8)
    with pytest.raises(ValueError, match="limit"):
        toads_refresh_relation_state(segmentation, memberships, boundary, limit)  # type: ignore[arg-type]


def test_relation_state_refresh_rejects_class_count_mismatch() -> None:
    segmentation = np.zeros((3, 3, 3), dtype=np.int8)
    memberships = np.ones((3, 3, 3, 2), dtype=np.float32)
    boundary = np.zeros((3, 3, 3, 3), dtype=np.int8)
    with pytest.raises(ValueError, match="same class count"):
        toads_refresh_relation_state(segmentation, memberships, boundary, 1)
