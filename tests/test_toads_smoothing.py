from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_math import toads_six_neighbor_smoothing_terms


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_smoothing"
SOURCE_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
SMOOTHING = np.float32(0.017)
PRIORS = np.asarray([0.5, 0.8, 1.1, 0.25], dtype=np.float32)


def _load_memberships() -> np.ndarray:
    image = nib.load(str(ROOT / "memberships.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj, dtype=np.float32)


def _expected_memberships() -> np.ndarray:
    x, y, z, class_index = np.indices((7, 6, 5, 4), dtype=np.int32)
    return (((7 * x + 11 * y + 13 * z + 5 * class_index) % 47) / 46.0).astype(
        np.float32
    )


def test_six_neighbor_smoothing_terms_match_java_reference() -> None:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["shape"] == [7, 6, 5, 4]
    assert reference["smoothing"] == 0.017
    np.testing.assert_array_equal(np.asarray(reference["smoothing_priors"], dtype=np.float32), PRIORS)
    memberships = _load_memberships()
    np.testing.assert_array_equal(memberships, _expected_memberships())

    num, energy = toads_six_neighbor_smoothing_terms(memberships, SMOOTHING, PRIORS)
    np.testing.assert_array_equal(num, np.asarray(reference["num"], dtype=np.float32).reshape(num.shape))
    np.testing.assert_array_equal(
        energy, np.asarray(reference["energy"], dtype=np.float32).reshape(energy.shape)
    )
    assert not np.any(num[[0, -1], :, :, :])
    assert not np.any(num[:, [0, -1], :, :])
    assert not np.any(num[:, :, [0, -1], :])


@pytest.mark.parametrize("smoothing", [0.0, -0.1])
def test_nonpositive_smoothing_matches_source_guard(smoothing: float) -> None:
    memberships = _expected_memberships()
    num, energy = toads_six_neighbor_smoothing_terms(memberships, smoothing, PRIORS)
    assert not np.any(num)
    assert not np.any(energy)


@pytest.mark.parametrize(
    "offset",
    [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)],
)
def test_each_face_neighbor_contributes_but_candidate_class_is_excluded(offset) -> None:
    memberships = np.full((5, 5, 5, 3), np.float32(0.2), dtype=np.float32)
    center = (2, 2, 2)
    neighbor = tuple(center[axis] + offset[axis] for axis in range(3))
    baseline_num, baseline_energy = toads_six_neighbor_smoothing_terms(
        memberships, np.float32(0.1), np.ones(3, dtype=np.float32)
    )

    changed = memberships.copy()
    changed[neighbor + (1,)] = np.float32(0.9)
    changed_num, changed_energy = toads_six_neighbor_smoothing_terms(
        changed, np.float32(0.1), np.ones(3, dtype=np.float32)
    )
    for candidate in (0, 2):
        assert changed_num[center + (candidate,)] > baseline_num[center + (candidate,)]
        assert changed_energy[center + (candidate,)] > baseline_energy[center + (candidate,)]
    assert changed_num[center + (1,)] == baseline_num[center + (1,)]
    assert changed_energy[center + (1,)] == baseline_energy[center + (1,)]


@pytest.mark.parametrize(
    ("memberships", "smoothing", "priors"),
    [
        (np.zeros((3, 3, 3), dtype=np.float32), 0.1, np.ones(1)),
        (np.zeros((3, 3, 3, 2), dtype=np.float32), 0.1, np.ones(1)),
        (np.zeros((3, 3, 3, 1), dtype=np.float32), np.nan, np.ones(1)),
        (np.full((3, 3, 3, 1), np.nan, dtype=np.float32), 0.1, np.ones(1)),
    ],
)
def test_six_neighbor_smoothing_rejects_invalid_inputs(memberships, smoothing, priors) -> None:
    with pytest.raises(ValueError):
        toads_six_neighbor_smoothing_terms(memberships, smoothing, priors)
