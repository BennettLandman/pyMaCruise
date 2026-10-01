from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_initial_relation_state import toads_initialize_relation_state


REFERENCE = (
    Path(__file__).parent
    / "data/synthetic/toads_initial_relation_state/java_reference.json"
)


def test_initial_relation_state_matches_selected_java_startup_sequence() -> None:
    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    assert reference["source_sha256"] == (
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    )
    assert reference["boundary_initialization"] == "zero"
    shape = tuple(reference["shape"])
    classes = reference["classes"]
    segmentation = np.asarray(reference["segmentation"], dtype=np.int8).reshape(shape)
    membership_bits = np.asarray(reference["membership_float_bits"], dtype=np.uint32)
    memberships = membership_bits.view(np.float32).reshape((*shape, classes))
    segmentation_before = segmentation.copy()
    memberships_before = memberships.copy()

    actual = toads_initialize_relation_state(
        segmentation, memberships, reference["limit"]
    )

    # Startup relation, boundary, and closest-label preparation reads the
    # caller's initialized state; it must not mutate either input in place.
    np.testing.assert_array_equal(segmentation, segmentation_before)
    np.testing.assert_array_equal(memberships, memberships_before)
    np.testing.assert_array_equal(actual.relation.relation_map.ravel(), reference["relation_map"])
    np.testing.assert_array_equal(
        actual.relation.relation_membership.ravel().astype(np.int8),
        reference["relation_membership"],
    )
    np.testing.assert_array_equal(actual.boundary_distances.ravel(), reference["boundary_distances"])
    np.testing.assert_array_equal(actual.closest_labels.ravel(), reference["closest_labels"])
    assert actual.relation.diffusion_iterations == reference["relation_diffusion_iterations"]
    assert actual.boundary_diffusion_passes == reference["boundary_diffusion_passes"]
    assert actual.closest_label_passes == 5 * reference["limit"] - 1


@pytest.mark.parametrize(
    ("segmentation", "memberships"),
    [
        (np.zeros((2, 3, 3), dtype=np.int8), np.zeros((2, 3, 3, 2), dtype=np.float32)),
        (np.zeros((3, 3, 3), dtype=np.float32), np.zeros((3, 3, 3, 2), dtype=np.float32)),
        (np.zeros((3, 3, 3), dtype=np.int8), np.zeros((2, 3, 3, 2), dtype=np.float32)),
    ],
)
def test_initial_relation_state_rejects_invalid_inputs(segmentation, memberships) -> None:
    with pytest.raises(ValueError):
        toads_initialize_relation_state(segmentation, memberships, 1)
