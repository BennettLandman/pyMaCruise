from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import iter_lesion_toads_volume_smoothing_updates


REFERENCE = (
    Path(__file__).parent
    / "data/synthetic/lesion_toads_volume_smoothing/java_reference.json"
)
UPDATE_FIELDS = (
    "class_membership_data_term",
    "lesion_membership_data_term",
    "black_hole_membership_data_term",
    "outmask_membership_data_term",
    "class_energy",
    "lesion_energy",
    "black_hole_energy",
    "other_class_neighbor_squares",
    "lesion_neighbor_squares",
    "black_hole_neighbor_squares",
    "current_class_neighbor_squares",
)


def _bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def _synthetic_inputs(shape: tuple[int, int, int], classes: int):
    x, y, z = np.indices(shape, dtype=np.float32)
    memberships = np.empty((*shape, classes), dtype=np.float32)
    class_terms = np.empty_like(memberships)
    lesion_terms = np.empty_like(memberships)
    black_hole_terms = np.empty_like(memberships)
    outmask_terms = np.empty_like(memberships)
    class_energies = np.empty_like(memberships)
    lesion_energies = np.empty_like(memberships)
    black_hole_energies = np.empty_like(memberships)
    for k in range(classes):
        memberships[..., k] = 1.0 + 0.5 * x + 0.25 * y + 0.125 * z + 0.0625 * k
        coordinate_sum = 1.0 + x + y + z + k
        class_terms[..., k] = 0.0625 * coordinate_sum
        lesion_terms[..., k] = 0.03125 * (coordinate_sum + 1.0)
        black_hole_terms[..., k] = 0.015625 * (coordinate_sum + 2.0)
        outmask_terms[..., k] = 0.125 * coordinate_sum
        class_energies[..., k] = 0.25 * coordinate_sum
        lesion_energies[..., k] = 0.125 * (coordinate_sum + 1.0)
        black_hole_energies[..., k] = 0.0625 * (coordinate_sum + 2.0)
    lesion = 0.25 + 0.125 * x + 0.0625 * y + 0.03125 * z
    black_hole = 0.5 + 0.0625 * x + 0.03125 * y + 0.015625 * z
    return (
        memberships,
        lesion.astype(np.float32),
        black_hole.astype(np.float32),
        class_terms,
        lesion_terms,
        black_hole_terms,
        outmask_terms,
        class_energies,
        lesion_energies,
        black_hole_energies,
    )


def test_streamed_smoothing_updates_match_selected_java_volume_harness() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    class_types = ("OBJECT", "LESION", "BLACKHOLE")
    priors = np.asarray([0.25, 0.5, 0.75], dtype=np.float32)
    inputs = _synthetic_inputs(tuple(reference["shape"]), len(class_types))
    before = tuple(value.copy() for value in inputs)
    result = list(
        iter_lesion_toads_volume_smoothing_updates(
            *inputs,
            class_types,
            priors,
            smoothing=float(np.asarray(reference["smoothing_bits"], dtype=np.uint32).view(np.float32)),
        )
    )

    assert [list(item.coordinate) for item in result] == [
        voxel["coordinate"] for voxel in reference["voxels"]
    ]
    for item, expected_voxel in zip(result, reference["voxels"], strict=True):
        assert len(item.class_updates) == len(expected_voxel["classes"]) == len(class_types)
        for update, expected in zip(item.class_updates, expected_voxel["classes"], strict=True):
            assert [_bits(getattr(update, name)) for name in UPDATE_FIELDS] == expected
    for initial, current in zip(before, inputs, strict=True):
        np.testing.assert_array_equal(initial, current)


def test_volume_smoothing_rejects_misaligned_class_term_maps() -> None:
    values = _synthetic_inputs((3, 3, 3), 2)
    wrong_terms = np.zeros((3, 3, 3, 1), dtype=np.float32)
    values = list(values)
    values[3] = wrong_terms
    with pytest.raises(ValueError, match="class_membership_data_terms must match"):
        list(
            iter_lesion_toads_volume_smoothing_updates(
                *values,
                ("OBJECT", "LESION"),
                np.asarray([0.25, 0.5]),
                smoothing=0.2,
            )
        )
