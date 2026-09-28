from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_math import toads_shape_weighted_membership_update


ROOT = Path(__file__).parent / "data" / "synthetic"
REFERENCE = ROOT / "toads_shape_weighted_membership" / "java_reference.json"
ZERO = np.float32(1e-10)
PRIOR_FACTOR = np.float32(0.75)


def _load(directory: str, name: str) -> np.ndarray:
    image = nib.load(str(ROOT / directory / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj)


def _inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    energy = _load("toads_direct_membership", "data_energy").astype(np.float32)
    initial = _load("toads_direct_membership", "initial_memberships").astype(np.float32)
    shapes = _load("toads_shape_relations", "shape_priors").astype(np.float32)
    labels = _load("toads_shape_relations", "closest_labels")
    return energy, initial, shapes, labels


def test_shape_weighted_membership_update_matches_java_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["shape"] == [7, 6, 5, 4]
    assert np.float32(reference["zero"]) == ZERO
    assert np.float32(reference["inf"]) == np.float32(1e10)
    assert np.float32(reference["prior_factor"]) == PRIOR_FACTOR
    energy, initial, shapes, labels = _inputs()
    updated = toads_shape_weighted_membership_update(
        energy, initial, shapes, labels, PRIOR_FACTOR
    )
    java = np.asarray(reference["updated"], dtype=np.float32).reshape(updated.shape)
    np.testing.assert_array_equal(updated, java)
    np.testing.assert_array_equal(updated[0], initial[0])
    np.testing.assert_array_equal(updated[:, -1], initial[:, -1])
    denominator = np.sum(updated[1:-1, 1:-1, 1:-1], axis=-1, dtype=np.float32)
    np.testing.assert_allclose(denominator, 1.0, rtol=0, atol=2e-7)


def test_zero_shape_prior_only_candidates_are_removed_before_normalization() -> None:
    energy = np.ones((3, 3, 3, 2), dtype=np.float32)
    initial = np.full_like(energy, 0.5)
    shapes = np.zeros_like(energy)
    labels = np.full((3, 3, 3, 4), -1, dtype=np.int16)
    labels[1, 1, 1, 0] = 0
    updated = toads_shape_weighted_membership_update(energy, initial, shapes, labels, 0.75)
    assert updated[1, 1, 1, 0] == 1.0
    assert updated[1, 1, 1, 1] == 0.0


def test_all_unsupported_classes_yield_zero_without_normalizing() -> None:
    energy = np.ones((3, 3, 3, 2), dtype=np.float32)
    initial = np.full_like(energy, 0.5)
    shapes = np.zeros_like(energy)
    labels = np.full((3, 3, 3, 4), -1, dtype=np.int16)
    updated = toads_shape_weighted_membership_update(energy, initial, shapes, labels, 0.75)
    assert np.all(updated[1, 1, 1] == 0.0)
    np.testing.assert_array_equal(updated[0], initial[0])


@pytest.mark.parametrize(
    ("energy", "initial", "shapes", "labels"),
    [
        (np.zeros((3, 3, 3), dtype=np.float32), np.zeros((3, 3, 3)), np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3, 4), dtype=np.int8)),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((2, 3, 3, 1)), np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3, 4), dtype=np.int8)),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((3, 3, 3, 1)), np.zeros((2, 3, 3, 1)), np.zeros((3, 3, 3, 4), dtype=np.int8)),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3, 1)), np.ones((3, 3, 3, 4), dtype=np.int8)),
    ],
)
def test_shape_weighted_membership_update_rejects_invalid_inputs(energy, initial, shapes, labels) -> None:
    with pytest.raises(ValueError):
        toads_shape_weighted_membership_update(energy, initial, shapes, labels, 0.75)
