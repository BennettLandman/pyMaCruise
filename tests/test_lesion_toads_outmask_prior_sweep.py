from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_membership import run_lesion_toads_outmask_smoothing_sweep


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_outmask_prior_sweep/java_reference.json"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
CLASS_TYPES = ("OBJECT", "OUTMASK", "OBJECT")
RELATIONS = ("direct", "prior", "direct")
SMOOTHING_PRIORS = np.asarray([0.25, 0.5, 0.75], dtype=np.float32)
ATLAS_PRIORS = np.asarray([0.5, 0.75, 1.25], dtype=np.float32)
SIMILARITIES = np.asarray([[0.0, 0.2, 0.7], [0.4, 0.0, 0.3], [0.6, 0.1, 0.0]], dtype=np.float32)


def _inputs(shape: tuple[int, int, int], classes: int):
    x, y, z = np.indices(shape, dtype=np.float32)
    memberships = np.empty((*shape, classes), dtype=np.float32)
    data = np.empty_like(memberships)
    shapes = np.empty_like(memberships)
    for k in range(classes):
        memberships[..., k] = 0.5 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        data[..., k] = 0.75 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        shapes[..., k] = 0.25 + 0.0625 * x + 0.03125 * y + 0.015625 * z + 0.0625 * k
    outmask_data = (0.375 + 0.03125 * x + 0.015625 * y + 0.0078125 * z).astype(np.float32)
    outlier = (0.125 * x + 0.0625 * y + 0.03125 * z).astype(np.float32)
    outlier[1, 1, 1] = 1.0
    if shape[0] > 3 and shape[1] > 2 and shape[2] > 4:
        outlier[3, 2, 4] = 0.0
    shapes[2, 2, 2, 1] = 0.0
    return memberships, outlier, data, outmask_data, shapes


def _run(inputs, *, coefficient: float, priors: bool):
    arguments = {}
    if priors:
        arguments = {
            "shape_priors": ATLAS_PRIORS,
            "class_similarities": SIMILARITIES,
            "prior_coefficient": coefficient,
        }
    return run_lesion_toads_outmask_smoothing_sweep(
        *inputs,
        SMOOTHING_PRIORS,
        CLASS_TYPES,
        RELATIONS,
        smoothing=0.2,
        **arguments,
    )


def test_outmask_prior_sweep_matches_selected_java_volume_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert int(np.asarray(np.float32(0.65)).view(np.uint32)) == reference["prior_coefficient_bits"]
    inputs = _inputs(tuple(reference["shape"]), reference["classes"])
    original_classes = inputs[0].copy()
    original_outlier = inputs[1].copy()
    result = _run(inputs, coefficient=float(np.float32(0.65)), priors=True)
    for expected in reference["voxels"]:
        point = tuple(expected["coordinate"])
        actual_classes = np.asarray(result.class_memberships[point], dtype=np.float32).view(np.uint32).tolist()
        assert actual_classes == expected["class_bits"]
        actual_outlier = int(np.asarray(result.outlier_membership[point], dtype=np.float32).view(np.uint32))
        assert actual_outlier == expected["outlier_bits"]
    border = np.ones(original_classes.shape[:3], dtype=bool)
    border[1:-1, 1:-1, 1:-1] = False
    np.testing.assert_array_equal(result.class_memberships[border], original_classes[border])
    np.testing.assert_array_equal(result.outlier_membership[border], original_outlier[border])


def test_zero_atlas_prior_coefficient_reproduces_outmask_smoothing_only() -> None:
    inputs = _inputs((5, 4, 6), 3)
    baseline = _run(inputs, coefficient=0.0, priors=False)
    with_zero_prior = _run(inputs, coefficient=0.0, priors=True)
    np.testing.assert_array_equal(with_zero_prior.class_memberships, baseline.class_memberships)
    np.testing.assert_array_equal(with_zero_prior.outlier_membership, baseline.outlier_membership)
