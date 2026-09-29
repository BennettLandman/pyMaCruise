from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import run_lesion_toads_lesion_outmask_smoothing_sweep


REFERENCE_ROOT = Path(__file__).parent / "data/synthetic/lesion_toads_lesion_outmask_distance_sweep"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
CLASS_TYPES = ("OBJECT", "LESION", "OUTMASK", "OBJECT")
RELATIONS = ("direct", "prior", "prior", "direct")
SMOOTHING_PRIORS = np.asarray([0.25, 0.5, 0.75, 0.375], dtype=np.float32)
ATLAS_PRIORS = np.asarray([0.5, 0.75, 1.25, 0.875], dtype=np.float32)
CLASS_SIMILARITIES = np.asarray(
    [[0.0, 0.2, 0.7, 0.3], [0.4, 0.0, 0.3, 0.6], [0.6, 0.1, 0.0, 0.2], [0.2, 0.5, 0.4, 0.0]],
    dtype=np.float32,
)
LESION_SIMILARITIES = np.asarray([0.25, 0.0, 0.5, 0.375], dtype=np.float32)


def _inputs(shape: tuple[int, int, int]):
    x, y, z = np.indices(shape, dtype=np.float32)
    memberships = np.empty((*shape, 4), dtype=np.float32)
    data = np.empty_like(memberships)
    shapes = np.empty_like(memberships)
    for k in range(4):
        memberships[..., k] = 0.5 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        data[..., k] = 0.75 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        shapes[..., k] = 0.25 + 0.0625 * x + 0.03125 * y + 0.015625 * z + 0.0625 * k
    lesion = (0.25 + 0.0625 * x + 0.03125 * y + 0.015625 * z).astype(np.float32)
    lesion_data = (0.625 + 0.03125 * x + 0.015625 * y + 0.0078125 * z).astype(np.float32)
    outmask_data = (0.375 + 0.03125 * x + 0.015625 * y + 0.0078125 * z).astype(np.float32)
    outlier = (0.125 * x + 0.0625 * y + 0.03125 * z).astype(np.float32)
    lesion[1, 1, 1] = 0.0
    outlier[1, 1, 1] = 1.0
    outlier[3, 2, 4] = 0.0
    shapes[2, 2, 2, 1] = 0.0
    ventricle = np.where(x == 1, 0, np.where(x == 2, 3, -1)).astype(np.int8)
    gray_matter = np.where(y == 1, 0, 4).astype(np.int8)
    brainstem = np.where(z == 2, 0, 2).astype(np.int8)
    return memberships, lesion, outlier, data, lesion_data, outmask_data, shapes, ventricle, gray_matter, brainstem


def _run(inputs, *, distances: bool, disabled_maps: bool = False, relations=RELATIONS):
    distance_arguments = {}
    if distances:
        distance_arguments = {
            "ventricle_distances": inputs[7],
            "gray_matter_distances": inputs[8],
            "brainstem_distances": inputs[9],
            "max_ventricle_distance": 5,
            "max_gray_matter_distance": 8,
            "max_brainstem_distance": 4,
            "gm_flag": True,
            "spread": 2.5,
        }
    elif disabled_maps:
        distance_arguments = {
            "ventricle_distances": inputs[7],
            "gray_matter_distances": inputs[8],
            "brainstem_distances": inputs[9],
            "gm_flag": True,
            "spread": 2.5,
        }
    return run_lesion_toads_lesion_outmask_smoothing_sweep(
        *inputs[:7],
        SMOOTHING_PRIORS,
        CLASS_TYPES,
        relations,
        smoothing=0.2,
        shape_priors=ATLAS_PRIORS,
        class_similarities=CLASS_SIMILARITIES,
        lesion_similarities=LESION_SIMILARITIES,
        prior_coefficient=0.65,
        **distance_arguments,
    )


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


@pytest.mark.parametrize(
    ("reference_name", "relations"),
    (("java_reference_prior.json", RELATIONS), ("java_reference_direct.json", ("direct", "direct", "prior", "direct"))),
)
def test_combined_distance_sweep_matches_selected_java_volume_reference(reference_name, relations) -> None:
    reference = json.loads((REFERENCE_ROOT / reference_name).read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["lesion_relation"] == relations[1]
    inputs = _inputs(tuple(reference["shape"]))
    original_classes, original_lesion, original_outlier = (array.copy() for array in inputs[:3])
    result = _run(inputs, distances=True, relations=relations)
    for expected in reference["voxels"]:
        point = tuple(expected["coordinate"])
        assert _bits(result.class_memberships[point]) == expected["class_bits"]
        assert int(np.asarray(result.lesion_membership[point], dtype=np.float32).view(np.uint32)) == expected["lesion_bits"]
        assert int(np.asarray(result.outlier_membership[point], dtype=np.float32).view(np.uint32)) == expected["outlier_bits"]
    border = np.ones(original_classes.shape[:3], dtype=bool)
    border[1:-1, 1:-1, 1:-1] = False
    np.testing.assert_array_equal(result.class_memberships[border], original_classes[border])
    np.testing.assert_array_equal(result.lesion_membership[border], original_lesion[border])
    np.testing.assert_array_equal(result.outlier_membership[border], original_outlier[border])


def test_disabling_combined_distance_maps_reproduces_p39zy() -> None:
    inputs = _inputs((5, 4, 6))
    without_maps = _run(inputs, distances=False)
    explicit_disabled_maps = _run(inputs, distances=False, disabled_maps=True)
    for actual, expected in zip(
        (without_maps.class_memberships, without_maps.lesion_membership, without_maps.outlier_membership),
        (explicit_disabled_maps.class_memberships, explicit_disabled_maps.lesion_membership, explicit_disabled_maps.outlier_membership),
        strict=True,
    ):
        np.testing.assert_array_equal(actual, expected)
