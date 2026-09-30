from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_objective import reduce_lesion_toads_objective


FIXTURE = Path(__file__).parent / "data/synthetic/lesion_objective_reduction"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
INPUT_SHA256 = "96d1ea8b4b739801067bb05fd04ebf7ee3f57ba63dd455387154717e2e21bb56"


def _bits(value: np.float32 | float) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def _reduce(volume: np.ndarray):
    return reduce_lesion_toads_objective(
        volume[..., :3],
        volume[..., 3:6],
        volume[..., 6],
        volume[..., 7],
        volume[..., 8],
        volume[..., 9],
        volume[..., 10],
    )


def test_objective_reduction_matches_selected_java_excerpt_on_ras_fixture() -> None:
    path = FIXTURE / "objective_inputs_ras.nii.gz"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == INPUT_SHA256
    image = nib.load(path)
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert image.shape == (5, 6, 7, 11)
    volume = np.asarray(image.dataobj, dtype=np.float32)
    assert np.isnan(volume[1, 1, 4, 6])
    assert np.isposinf(volume[3, 4, 3, 5])
    assert np.isnan(volume[3, 4, 5, 4])

    reference = json.loads((FIXTURE / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["shape"] == list(image.shape[:3])
    assert reference["classes"] == 3
    assert len(reference["voxels"]) == (5 - 2) * (6 - 2) * (7 - 2)
    result = _reduce(volume)
    for expected in reference["voxels"]:
        coordinate = tuple(expected["coordinate"])
        assert _bits(result.running_total[coordinate]) == expected["running_total_bits"], coordinate
    assert _bits(result.total) == reference["total_bits"]
    boundary = np.ones(image.shape[:3], dtype=bool)
    boundary[1:-1, 1:-1, 1:-1] = False
    np.testing.assert_array_equal(result.running_total[boundary], np.zeros(np.count_nonzero(boundary), dtype=np.float32))


def test_objective_reduction_preserves_java_class_loop_order() -> None:
    memberships = np.zeros((3, 3, 3, 3), dtype=np.float32)
    energies = np.zeros_like(memberships)
    denominator = np.zeros((3, 3, 3), dtype=np.float32)
    lesion = np.zeros_like(denominator)
    black_hole = np.zeros_like(denominator)
    lesion_energy = np.zeros_like(denominator)
    black_hole_energy = np.zeros_like(denominator)
    memberships[1, 1, 1] = [0.5, 0.25, 0.25]
    energies[1, 1, 1] = [2**26, 16.0, 16.0]
    denominator[1, 1, 1] = 1.0
    forward = reduce_lesion_toads_objective(
        memberships, energies, denominator, lesion, black_hole, lesion_energy, black_hole_energy
    )
    order = [1, 2, 0]
    reordered = reduce_lesion_toads_objective(
        memberships[..., order], energies[..., order], denominator, lesion, black_hole, lesion_energy, black_hole_energy
    )
    assert forward.total == np.float32(2**24)
    assert reordered.total == np.float32(2**24 + 2)
    assert forward.total != reordered.total


def test_auxiliary_energy_is_added_even_when_denominator_is_zero() -> None:
    shape = (3, 3, 3)
    memberships = np.ones((*shape, 2), dtype=np.float32)
    class_energies = np.full_like(memberships, 100.0)
    denominator = np.zeros(shape, dtype=np.float32)
    lesion = np.zeros(shape, dtype=np.float32)
    black_hole = np.zeros(shape, dtype=np.float32)
    lesion_energy = np.zeros(shape, dtype=np.float32)
    black_hole_energy = np.zeros(shape, dtype=np.float32)
    lesion[1, 1, 1] = 0.5
    black_hole[1, 1, 1] = 0.25
    lesion_energy[1, 1, 1] = 2.0
    black_hole_energy[1, 1, 1] = 4.0
    result = reduce_lesion_toads_objective(
        memberships,
        class_energies,
        denominator,
        lesion,
        black_hole,
        lesion_energy,
        black_hole_energy,
    )
    assert result.total == np.float32(0.75)
