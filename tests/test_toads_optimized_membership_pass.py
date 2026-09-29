from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_math import (
    toads_optimized_centroid_similarity,
    toads_optimized_class_standard_deviation,
    toads_optimized_image_variance,
    toads_optimized_membership_pass,
    toads_optimized_tissue_data_term,
    toads_psf_ridge_bias,
)


ROOT = Path(__file__).parent / "data" / "synthetic"
REFERENCE_ROOT = ROOT / "toads_optimized_membership_pass"
SEGMENTATION_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
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


def _reference(field: bool, psf: bool) -> dict:
    suffix = f"field_{'on' if field else 'off'}_psf_{'on' if psf else 'off'}"
    return json.loads((REFERENCE_ROOT / f"java_reference_{suffix}.json").read_text())


@pytest.mark.parametrize("use_field", [False, True])
@pytest.mark.parametrize("use_psf", [False, True])
def test_composed_optimized_tissue_membership_pass_matches_java(
    use_field: bool, use_psf: bool
) -> None:
    reference = _reference(use_field, use_psf)
    assert reference["segmentation_source_sha256"] == SEGMENTATION_SHA256
    images, field, memberships, segmentation, shapes, labels = _inputs()
    field_values = field if use_field else None
    ridge = np.float32(0.2 if use_psf else 0.0)

    image_weights, image_stddev = toads_optimized_image_variance(
        images, memberships, CENTROIDS, intensity_field=field_values
    )
    class_stddev = toads_optimized_class_standard_deviation(
        images, memberships, CENTROIDS, segmentation, intensity_field=field_values
    )
    similarity = toads_optimized_centroid_similarity(
        CENTROIDS, class_stddev, image_weights, prior_scale=0.05
    )
    bias = toads_psf_ridge_bias(
        images, segmentation, ridge, intensity_field=field_values
    )
    data_energy = toads_optimized_tissue_data_term(
        images,
        CENTROIDS,
        EXPONENTS,
        MODALITY_WEIGHTS,
        CLASS_VARIANCE_WEIGHTS,
        image_weights,
        intensity_field=field_values,
        psf_bias=bias,
    )
    np.testing.assert_array_equal(image_weights, np.asarray(reference["image_variance_weights"], dtype=np.float32))
    np.testing.assert_array_equal(image_stddev, np.asarray(reference["image_stddev"], dtype=np.float32))
    np.testing.assert_array_equal(class_stddev, np.asarray(reference["class_standard_deviations"], dtype=np.float32).reshape(3, 3))
    np.testing.assert_array_equal(similarity, np.asarray(reference["similarity"], dtype=np.float32).reshape(3, 3))
    np.testing.assert_array_equal(bias, np.asarray(reference["bias"], dtype=np.float32).reshape(bias.shape))
    np.testing.assert_allclose(
        data_energy,
        np.asarray(reference["data_energy"], dtype=np.float32).reshape(data_energy.shape),
        rtol=0,
        atol=1e-7,
    )

    updated, objective = toads_optimized_membership_pass(
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
        intensity_field=field_values,
        psf_ridge=ridge,
        prior_scale=0.05,
    )
    np.testing.assert_allclose(
        updated,
        np.asarray(reference["updated"], dtype=np.float32).reshape(updated.shape),
        rtol=0,
        atol=2e-6,
    )
    assert np.isclose(objective, np.float32(reference["objective"]), rtol=0, atol=2e-6)
    np.testing.assert_array_equal(updated[0], memberships[0])
    np.testing.assert_array_equal(updated[:, -1], memberships[:, -1])


def test_composed_pass_rejects_missing_closest_label_slots() -> None:
    images, field, memberships, segmentation, shapes, labels = _inputs()
    with pytest.raises(ValueError, match="four label slots"):
        toads_optimized_membership_pass(
            images,
            memberships,
            CENTROIDS,
            EXPONENTS,
            MODALITY_WEIGHTS,
            CLASS_VARIANCE_WEIGHTS,
            segmentation,
            shapes,
            labels[..., :3],
            ATLAS_SHAPE_PRIORS,
            SMOOTHING_PRIORS,
            smoothing=0.017,
            prior_factor=0.75,
            prior_coefficient=0.15,
        )
