from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_math import toads_direct_membership_update


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_direct_membership"
SOURCE_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
ZERO = np.float32(1e-10)


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj, dtype=np.float32)


def _expected_inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    shape = (7, 6, 5, 4)
    x, y, z, class_index = np.indices(shape, dtype=np.int32)
    energy = (((13 * x + 7 * y + 3 * z + 11 * class_index) % 29) + 1).astype(np.float32) / np.float32(13.0)
    objective = (((5 * x + 3 * y + 7 * z + 2 * class_index) % 19) + 1).astype(np.float32) / np.float32(17.0)
    initial = (((5 * x + 7 * y + 11 * z + 3 * class_index) % 31) / 30.0).astype(np.float32)
    energy[2, 2, 2, 0] = 0.0
    energy[3, 2, 2, 1] = ZERO
    energy[4, 3, 2, 2] = np.nextafter(ZERO, np.float32(0.0))
    energy[2, 3, 2, 3] = np.nextafter(ZERO, np.float32(1.0))
    return energy, objective, initial


def test_direct_membership_update_matches_java_reference() -> None:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["shape"] == [7, 6, 5, 4]
    assert np.float32(reference["zero"]) == ZERO
    assert np.float32(reference["inf"]) == np.float32(1e10)
    energy, objective, initial = (_load("data_energy"), _load("objective_energy"), _load("initial_memberships"))
    expected = _expected_inputs()
    for actual, formula in zip((energy, objective, initial), expected, strict=True):
        np.testing.assert_array_equal(actual, formula)

    updated, contributions = toads_direct_membership_update(
        energy, initial, objective_energy=objective
    )
    java_updated = np.asarray(reference["updated"], dtype=np.float32).reshape(updated.shape)
    java_contributions = np.asarray(reference["objective_contributions"], dtype=np.float32).reshape(
        contributions.shape
    )
    np.testing.assert_array_equal(updated, java_updated)
    np.testing.assert_array_equal(contributions, java_contributions)
    np.testing.assert_array_equal(updated[0], initial[0])
    np.testing.assert_array_equal(updated[:, -1], initial[:, -1])
    np.testing.assert_array_equal(contributions[[0, -1]], 0.0)

    # ZERO, below-ZERO, exact ZERO, and above-ZERO exercise source `>` logic.
    for coordinate, winner in (
        ((2, 2, 2), 0),
        ((3, 2, 2), 1),
        ((4, 3, 2), 2),
        ((2, 3, 2), 3),
    ):
        assert int(np.argmax(updated[coordinate])) == winner


@pytest.mark.parametrize(
    ("energy", "initial", "objective", "zero_threshold", "infinite_reciprocal"),
    [
        (np.zeros((3, 3, 3), dtype=np.float32), np.zeros((3, 3, 3)), None, 1e-10, 1e10),
        (np.zeros((3, 3, 3, 2), dtype=np.float32), np.zeros((2, 3, 3, 2)), None, 1e-10, 1e10),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3)), 1e-10, 1e10),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((3, 3, 3, 1)), None, -1.0, 1e10),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((3, 3, 3, 1)), None, 1e-10, 0.0),
    ],
)
def test_direct_membership_update_rejects_invalid_inputs(
    energy, initial, objective, zero_threshold, infinite_reciprocal
) -> None:
    with pytest.raises(ValueError):
        toads_direct_membership_update(
            energy,
            initial,
            objective_energy=objective,
            zero_threshold=zero_threshold,
            infinite_reciprocal=infinite_reciprocal,
        )
