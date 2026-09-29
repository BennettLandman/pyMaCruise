from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_math import toads_inverse_shape_prior_energy


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_inverse_shape_prior"
REFERENCE = ROOT / "java_reference.json"
SEGMENTATION_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def _inputs() -> tuple[dict, np.ndarray]:
    reference = json.loads(REFERENCE.read_text())
    image = nib.load(str(ROOT / "shape_values.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return reference, np.asarray(image.dataobj, dtype=np.float32)


def test_inverse_shape_prior_energy_matches_java_formula_reference() -> None:
    reference, shapes = _inputs()
    assert reference["segmentation_source_sha256"] == SEGMENTATION_SHA256
    assert list(shapes.shape) == reference["shape"]
    inverse, objective = toads_inverse_shape_prior_energy(
        shapes,
        reference["shape_priors"],
        reference["similarity"],
        reference["prior_coefficient"],
    )
    np.testing.assert_array_equal(
        inverse, np.asarray(reference["inverse_energy"], dtype=np.float32).reshape(shapes.shape)
    )
    np.testing.assert_array_equal(
        objective,
        np.asarray(reference["objective_energy"], dtype=np.float32).reshape(shapes.shape),
    )
    assert np.all(inverse[[0, -1], ...] == 0.0)
    assert np.all(objective[:, [0, -1], ...] == 0.0)
    assert np.all(inverse[:, :, [0, -1], :] == 0.0)


def test_nonpositive_prior_coefficient_adds_no_terms() -> None:
    reference, shapes = _inputs()
    inverse, objective = toads_inverse_shape_prior_energy(
        shapes,
        reference["shape_priors"],
        reference["similarity"],
        0.0,
    )
    assert not np.any(inverse)
    assert not np.any(objective)


def test_candidate_shape_is_excluded_from_its_own_inverse_prior() -> None:
    reference, shapes = _inputs()
    original, _ = toads_inverse_shape_prior_energy(
        shapes, reference["shape_priors"], reference["similarity"], reference["prior_coefficient"]
    )
    changed = shapes.copy()
    changed[2, 2, 2, 1] = np.float32(1.0) - changed[2, 2, 2, 1]
    updated, _ = toads_inverse_shape_prior_energy(
        changed, reference["shape_priors"], reference["similarity"], reference["prior_coefficient"]
    )
    assert updated[2, 2, 2, 1] == original[2, 2, 2, 1]


@pytest.mark.parametrize(
    ("shapes", "priors", "similarity"),
    [
        (np.zeros((3, 3, 3)), np.zeros(1), np.zeros((1, 1))),
        (np.zeros((3, 3, 3, 1)), np.zeros(2), np.zeros((1, 1))),
        (np.zeros((3, 3, 3, 1)), np.zeros(1), np.zeros((2, 2))),
    ],
)
def test_inverse_shape_prior_rejects_invalid_shapes(shapes, priors, similarity) -> None:
    with pytest.raises(ValueError):
        toads_inverse_shape_prior_energy(shapes, priors, similarity, 0.1)
