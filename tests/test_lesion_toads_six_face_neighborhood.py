from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import iter_lesion_toads_six_face_neighborhoods


REFERENCE = (
    Path(__file__).parent
    / "data/synthetic/lesion_toads_six_face_neighborhood/java_reference.json"
)


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def _synthetic_memberships(shape: tuple[int, int, int], classes: int):
    nx, ny, nz = shape
    x, y, z = np.indices(shape, dtype=np.float32)
    class_values = np.empty((*shape, classes), dtype=np.float32)
    for k in range(classes):
        class_values[..., k] = 1.0 + 0.5 * x + 0.25 * y + 0.125 * z + 0.0625 * k
    lesion = 0.25 + 0.125 * x + 0.0625 * y + 0.03125 * z
    black_hole = 0.5 + 0.0625 * x + 0.03125 * y + 0.015625 * z
    return class_values, lesion.astype(np.float32), black_hole.astype(np.float32)


def test_six_face_stencils_match_selected_java_loop_order_and_values() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert reference["face_order"] == ["X+", "X-", "Y+", "Y-", "Z+", "Z-"]
    classes, lesion, black_hole = _synthetic_memberships(tuple(reference["shape"]), reference["classes"])
    originals = (classes.copy(), lesion.copy(), black_hole.copy())

    neighborhoods = list(iter_lesion_toads_six_face_neighborhoods(classes, lesion, black_hole))
    assert len(neighborhoods) == len(reference["voxels"]) == (
        (classes.shape[0] - 2) * (classes.shape[1] - 2) * (classes.shape[2] - 2)
    )
    assert [list(item.coordinate) for item in neighborhoods] == [
        voxel["coordinate"] for voxel in reference["voxels"]
    ]
    for item, expected in zip(neighborhoods, reference["voxels"], strict=True):
        assert _bits(item.class_memberships) == expected["class_bits"]
        assert _bits(item.lesion_memberships) == expected["lesion_bits"]
        assert _bits(item.black_hole_memberships) == expected["black_hole_bits"]
    for before, after in zip(originals, (classes, lesion, black_hole), strict=True):
        np.testing.assert_array_equal(after, before)


def test_six_face_stencils_reject_invalid_volume_shapes() -> None:
    with pytest.raises(ValueError, match="spatial size >= 3"):
        list(iter_lesion_toads_six_face_neighborhoods(np.zeros((2, 3, 3, 1)), np.zeros((2, 3, 3)), np.zeros((2, 3, 3))))
    with pytest.raises(ValueError, match="match the class spatial shape"):
        list(iter_lesion_toads_six_face_neighborhoods(np.zeros((3, 3, 3, 1)), np.zeros((3, 3, 2)), np.zeros((3, 3, 3))))

