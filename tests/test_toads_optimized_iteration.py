from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_math import toads_optimized_tissue_iteration


ROOT = Path(__file__).parent / "data" / "synthetic"
REFERENCE_ROOT = ROOT / "toads_optimized_membership_pass"
CENTROIDS = np.asarray(
    [[0.64, 0.62, 0.81], [0.35, 0.55, 0.72], [0.18, 0.70, 0.88]],
    dtype=np.float32,
)
EXPONENTS = np.asarray([[0.8, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], dtype=np.float32)
MODALITY_WEIGHTS = np.asarray([1.0, 0.75, 0.5], dtype=np.float32)
CLASS_VARIANCE_WEIGHTS = np.ones((3, 3), dtype=np.float32)
ATLAS_SHAPE_PRIORS = np.asarray([0.2, 0.35, 0.45], dtype=np.float32)
SMOOTHING_PRIORS = np.asarray([0.5, 0.8, 1.1], dtype=np.float32)


def _load(directory: str, name: str, *, integer: bool = False) -> np.ndarray:
    image = nib.load(str(ROOT / directory / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj, dtype=np.int16 if integer else np.float32)


def _inputs() -> tuple[np.ndarray, ...]:
    images = _load("toads_optimized_likelihood", "images")
    field = _load("toads_optimized_likelihood", "field")
    memberships = _load("toads_direct_membership", "initial_memberships")[..., :3]
    segmentation = _load("toads_optimized_class_variance", "segmentation", integer=True)
    shapes = _load("toads_shape_relations", "shape_priors")[..., :3]
    labels = _load("toads_shape_relations", "closest_labels", integer=True)
    labels = np.where(labels < 3, labels, -1).astype(np.int16)
    return images, field, memberships, segmentation, shapes, labels


def _reference(use_field: bool, use_psf: bool) -> dict:
    suffix = f"field_{'on' if use_field else 'off'}_psf_{'on' if use_psf else 'off'}"
    return json.loads((REFERENCE_ROOT / f"java_iteration_{suffix}.json").read_text())


@pytest.mark.parametrize("use_field", [False, True])
@pytest.mark.parametrize("use_psf", [False, True])
def test_optimized_tissue_iteration_matches_java_loop_order(use_field: bool, use_psf: bool) -> None:
    reference = _reference(use_field, use_psf)
    assert reference["segmentation_source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert reference["power_table_source_sha256"] == "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
    images, field, memberships, segmentation, shapes, labels = _inputs()
    result, objective, centroids = toads_optimized_tissue_iteration(
        images,
        memberships,
        CENTROIDS,
        EXPONENTS,
        MODALITY_WEIGHTS,
        CLASS_VARIANCE_WEIGHTS,
        segmentation,
        shapes,
        labels,
        ATLAS_SHAPE_PRIORS,
        SMOOTHING_PRIORS,
        smoothing=0.017,
        prior_factor=0.75,
        prior_coefficient=0.15,
        centroid_smoothness=0.1,
        intensity_field=field if use_field else None,
        psf_ridge=0.2 if use_psf else 0.0,
    )

    np.testing.assert_array_equal(centroids, np.asarray(reference["centroids"], dtype=np.float32).reshape(3, 3))
    np.testing.assert_allclose(
        result,
        np.asarray(reference["updated"], dtype=np.float32).reshape(result.shape),
        rtol=0,
        atol=2e-6,
    )
    assert np.isclose(objective, np.float32(reference["objective"]), rtol=0, atol=2e-6)
    np.testing.assert_array_equal(result[0], memberships[0])
    np.testing.assert_array_equal(result[:, -1], memberships[:, -1])
