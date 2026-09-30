from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_relations import (
    lesion_build_relation_list,
    lesion_compute_relations,
)


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_relations/java_reference.json"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"


def test_lesion_relation_list_and_diffusion_match_selected_java_excerpts() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert len(reference["cases"]) == 6
    for case in reference["cases"]:
        shape = tuple(case["shape"])
        labels = np.asarray(case["segmentation"], dtype=np.int8).reshape(shape)
        before = labels.copy()

        built = lesion_build_relation_list(labels, case["classes"])
        expected_initial = np.asarray(case["initial_map"], dtype=np.int16).reshape(shape)
        np.testing.assert_array_equal(built.relation_map, expected_initial)
        np.testing.assert_array_equal(
            built.relation_membership.ravel(),
            np.asarray(case["relations"], dtype=np.bool_),
        )
        assert built.relation_count == case["relation_count"]
        assert built.diffusion_iterations == 0

        result = lesion_compute_relations(labels, case["classes"], case["limit"])
        expected = np.asarray(case["relation_map"], dtype=np.int16).reshape(shape)
        np.testing.assert_array_equal(result.relation_map, expected)
        np.testing.assert_array_equal(result.relation_membership, built.relation_membership)
        assert result.diffusion_iterations == case["iterations"]
        # Every border remains the constructor's original class ID, untouched
        # by both buildRelationList and computeRelations.
        for axis in range(3):
            for border in (0, -1):
                index = [slice(None)] * 3
                index[axis] = border
                np.testing.assert_array_equal(result.relation_map[tuple(index)], labels[tuple(index)])
        np.testing.assert_array_equal(labels, before)


def test_lesion_relation_diffusion_uses_first_source_face_and_synchronous_passes() -> None:
    reference = json.loads(REFERENCE.read_text())
    multi = {case["name"]: case for case in reference["cases"]}
    once = multi["asymmetric-planar-one-pass"]
    twice = multi["asymmetric-planar-two-passes"]
    first_map = np.asarray(once["relation_map"], dtype=np.int16).reshape(once["shape"])
    second_map = np.asarray(twice["relation_map"], dtype=np.int16).reshape(twice["shape"])
    initial = np.asarray(twice["initial_map"], dtype=np.int16).reshape(twice["shape"])
    classes = twice["classes"]
    # The two-pass map differs at newly reached class-only voxels, which checks
    # that a pass reads the preceding map rather than updating in place.
    assert np.count_nonzero((first_map < classes) & (second_map >= classes)) > 0
    assert np.count_nonzero(initial != first_map) > 0


@pytest.mark.parametrize(
    ("labels", "classes", "limit"),
    [
        (np.zeros((2, 4, 4), dtype=np.int8), 1, 1),
        (np.zeros((3, 3, 3), dtype=float), 1, 1),
        (np.full((3, 3, 3), -1, dtype=np.int8), 1, 1),
        (np.full((3, 3, 3), 1, dtype=np.int8), 1, 1),
        (np.zeros((3, 3, 3), dtype=np.int8), 0, 1),
        (np.zeros((3, 3, 3), dtype=np.int8), 1, -1),
    ],
)
def test_lesion_relations_reject_invalid_input(labels, classes, limit) -> None:
    with pytest.raises(ValueError):
        lesion_compute_relations(labels, classes, limit)


def test_lesion_relation_list_rejects_class_labels_outside_java_byte_range() -> None:
    labels = np.zeros((3, 3, 3), dtype=np.int8)
    with pytest.raises(ValueError, match="signed-byte"):
        lesion_build_relation_list(labels, 129)
