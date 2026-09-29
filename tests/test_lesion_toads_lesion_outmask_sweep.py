from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import run_lesion_toads_lesion_outmask_smoothing_sweep


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_lesion_outmask_sweep/java_reference.json"
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


def _inputs(shape: tuple[int, int, int], classes: int):
    x, y, z = np.indices(shape, dtype=np.float32)
    memberships = np.empty((*shape, classes), dtype=np.float32)
    data = np.empty_like(memberships)
    shapes = np.empty_like(memberships)
    for k in range(classes):
        memberships[..., k] = 0.5 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        data[..., k] = 0.75 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        shapes[..., k] = 0.25 + 0.0625 * x + 0.03125 * y + 0.015625 * z + 0.0625 * k
    lesion = (0.25 + 0.0625 * x + 0.03125 * y + 0.015625 * z).astype(np.float32)
    lesion_data = (0.625 + 0.03125 * x + 0.015625 * y + 0.0078125 * z).astype(np.float32)
    outmask_data = (0.375 + 0.03125 * x + 0.015625 * y + 0.0078125 * z).astype(np.float32)
    outlier = (0.125 * x + 0.0625 * y + 0.03125 * z).astype(np.float32)
    lesion[1, 1, 1] = 0.0
    outlier[1, 1, 1] = 1.0
    if shape[0] > 3 and shape[1] > 2 and shape[2] > 4:
        outlier[3, 2, 4] = 0.0
    shapes[2, 2, 2, 1] = 0.0
    return memberships, lesion, outlier, data, lesion_data, outmask_data, shapes


def _run(inputs, *, include_priors: bool, coefficient: float = 0.65):
    prior_arguments = {}
    if include_priors:
        prior_arguments = {
            "shape_priors": ATLAS_PRIORS,
            "class_similarities": CLASS_SIMILARITIES,
            "lesion_similarities": LESION_SIMILARITIES,
            "prior_coefficient": coefficient,
        }
    return run_lesion_toads_lesion_outmask_smoothing_sweep(
        *inputs,
        SMOOTHING_PRIORS,
        CLASS_TYPES,
        RELATIONS,
        smoothing=0.2,
        **prior_arguments,
    )


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def test_lesion_outmask_sweep_matches_selected_java_volume_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    inputs = _inputs(tuple(reference["shape"]), reference["classes"])
    original_classes, original_lesion, original_outlier = (array.copy() for array in inputs[:3])
    result = _run(inputs, include_priors=True)
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
    for actual, original in zip(inputs[:3], (original_classes, original_lesion, original_outlier), strict=True):
        np.testing.assert_array_equal(actual, original)


def test_zero_prior_combined_sweep_matches_prior_enabled_zero_control() -> None:
    inputs = _inputs((5, 4, 6), 4)
    no_atlas_prior = _run(inputs, include_priors=False)
    zero_atlas_prior = _run(inputs, include_priors=True, coefficient=0.0)
    np.testing.assert_array_equal(no_atlas_prior.class_memberships, zero_atlas_prior.class_memberships)
    np.testing.assert_array_equal(no_atlas_prior.lesion_membership, zero_atlas_prior.lesion_membership)
    np.testing.assert_array_equal(no_atlas_prior.outlier_membership, zero_atlas_prior.outlier_membership)


def test_lesion_outmask_sweep_rejects_duplicate_auxiliary_types() -> None:
    inputs = _inputs((3, 3, 3), 4)
    with pytest.raises(ValueError, match="exactly one LESION and one OUTMASK"):
        run_lesion_toads_lesion_outmask_smoothing_sweep(
            *inputs,
            SMOOTHING_PRIORS,
            ("OBJECT", "LESION", "OUTMASK", "OUTMASK"),
            RELATIONS,
            smoothing=0.2,
        )
