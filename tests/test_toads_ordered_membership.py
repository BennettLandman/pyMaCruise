from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_math import (
    toads_ordered_membership_sweep,
    toads_shape_weighted_membership_update,
    toads_six_neighbor_smoothing_terms,
)


ROOT = Path(__file__).parent / "data" / "synthetic"
REFERENCE = ROOT / "toads_ordered_membership" / "java_reference.json"
SMOOTHING = np.float32(0.017)
PRIOR_FACTOR = np.float32(0.75)
SMOOTHING_PRIORS = np.asarray([0.5, 0.8, 1.1, 0.25], dtype=np.float32)


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


def test_ordered_membership_sweep_matches_java_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["shape"] == [7, 6, 5, 4]
    assert np.float32(reference["smoothing"]) == SMOOTHING
    assert np.float32(reference["prior_factor"]) == PRIOR_FACTOR
    np.testing.assert_array_equal(
        np.asarray(reference["smoothing_priors"], dtype=np.float32), SMOOTHING_PRIORS
    )
    energy, initial, shapes, labels = _inputs()
    updated, objective = toads_ordered_membership_sweep(
        energy,
        initial,
        shapes,
        labels,
        SMOOTHING_PRIORS,
        SMOOTHING,
        PRIOR_FACTOR,
    )
    java = np.asarray(reference["updated"], dtype=np.float32).reshape(updated.shape)
    np.testing.assert_array_equal(updated, java)
    assert objective == np.float32(reference["objective"])
    np.testing.assert_array_equal(updated[0], initial[0])
    np.testing.assert_array_equal(updated[:, -1], initial[:, -1])
    np.testing.assert_allclose(
        np.sum(updated[1:-1, 1:-1, 1:-1], axis=-1, dtype=np.float32),
        1.0,
        rtol=0,
        atol=2e-7,
    )

    # Frozen-input smoothing makes a different update than the source's in-place scan.
    smooth_num, _ = toads_six_neighbor_smoothing_terms(initial, SMOOTHING, SMOOTHING_PRIORS)
    synchronous = toads_shape_weighted_membership_update(
        np.asarray(energy + smooth_num, dtype=np.float32),
        initial,
        shapes,
        labels,
        PRIOR_FACTOR,
    )
    assert not np.array_equal(updated, synchronous)


@pytest.mark.parametrize(
    ("energy", "initial", "shapes", "labels", "smoothing_priors"),
    [
        (np.zeros((3, 3, 3), dtype=np.float32), np.zeros((3, 3, 3)), np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3, 4), dtype=np.int8), np.ones(1)),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((2, 3, 3, 1)), np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3, 4), dtype=np.int8), np.ones(1)),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3, 4), dtype=np.int8), np.ones(2)),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 3, 1)), np.ones((3, 3, 3, 4), dtype=np.int8), np.ones(1)),
    ],
)
def test_ordered_membership_sweep_rejects_invalid_inputs(
    energy, initial, shapes, labels, smoothing_priors
) -> None:
    with pytest.raises(ValueError):
        toads_ordered_membership_sweep(
            energy, initial, shapes, labels, smoothing_priors, 0.017, 0.75
        )
