from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import run_lesion_toads_blackhole_smoothing_sweep


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_blackhole_sweep/java_reference.json"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
CLASS_TYPES = ("OBJECT", "BLACKHOLE", "OBJECT")
RELATIONS = ("direct", "prior", "direct")
SMOOTHING_PRIORS = np.asarray([0.25, 0.5, 0.75], dtype=np.float32)


def _inputs(shape: tuple[int, int, int], classes: int):
    x, y, z = np.indices(shape, dtype=np.float32)
    memberships = np.empty((*shape, classes), dtype=np.float32)
    class_data = np.empty_like(memberships)
    shapes = np.empty_like(memberships)
    for k in range(classes):
        memberships[..., k] = 0.5 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        class_data[..., k] = 0.75 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        shapes[..., k] = 0.25 + 0.0625 * x + 0.03125 * y + 0.015625 * z + 0.0625 * k
    lesion = (0.25 + 0.125 * x + 0.0625 * y + 0.03125 * z).astype(np.float32)
    black_hole = (0.5 + 0.0625 * x + 0.03125 * y + 0.015625 * z).astype(np.float32)
    lesion_data = (0.75 + 0.0625 * x + 0.03125 * y + 0.015625 * z).astype(np.float32)
    black_hole_data = (0.5 + 0.03125 * x + 0.015625 * y + 0.0078125 * z).astype(np.float32)
    shapes[2, 2, 2, 1] = 0.0
    return memberships, lesion, black_hole, class_data, lesion_data, black_hole_data, shapes


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def test_blackhole_sweep_matches_selected_java_volume_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    inputs = _inputs(tuple(reference["shape"]), reference["classes"])
    before = tuple(value.copy() for value in inputs[:3])
    result = run_lesion_toads_blackhole_smoothing_sweep(
        *inputs,
        SMOOTHING_PRIORS,
        CLASS_TYPES,
        RELATIONS,
        smoothing=float(np.asarray(reference["smoothing_bits"], dtype=np.uint32).view(np.float32)),
    )
    for expected in reference["voxels"]:
        point = tuple(expected["coordinate"])
        assert _bits(result.class_memberships[point]) == expected["class_bits"]
        assert int(np.asarray(result.lesion_membership[point], dtype=np.float32).view(np.uint32)) == expected["lesion_bits"]
        assert int(np.asarray(result.black_hole_membership[point], dtype=np.float32).view(np.uint32)) == expected["black_hole_bits"]
    border = np.ones(inputs[0].shape[:3], dtype=bool)
    border[1:-1, 1:-1, 1:-1] = False
    np.testing.assert_array_equal(result.class_memberships[border], before[0][border])
    np.testing.assert_array_equal(result.lesion_membership[border], before[1][border])
    np.testing.assert_array_equal(result.black_hole_membership[border], before[2][border])
    np.testing.assert_array_equal(inputs[0], before[0])
    np.testing.assert_array_equal(inputs[1], before[1])
    np.testing.assert_array_equal(inputs[2], before[2])


def test_blackhole_sweep_rejects_multiple_blackhole_classes() -> None:
    inputs = _inputs((3, 3, 3), 3)
    with pytest.raises(ValueError, match="exactly one BLACKHOLE class"):
        run_lesion_toads_blackhole_smoothing_sweep(
            *inputs,
            SMOOTHING_PRIORS,
            ("OBJECT", "BLACKHOLE", "BLACKHOLE"),
            RELATIONS,
            smoothing=0.2,
        )
