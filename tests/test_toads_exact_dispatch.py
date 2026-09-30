from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_exact_dispatch import toads_exact_membership_dispatch
from pymacruise.toads_math import (
    toads_object_mask_membership_pass,
    toads_object_mask_outlier_membership_pass,
)


ROOT = Path(__file__).parent / "data" / "synthetic"
REFERENCE = ROOT / "toads_object_mask_outlier_membership_pass" / "java_reference.json"
FOLDED = ROOT / "toads_exact_dispatch"
SEGMENTATION_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def _float(bits: int | list[int], shape: tuple[int, ...] | None = None) -> np.ndarray:
    array = np.asarray(bits, dtype=np.int32).view(np.float32)
    return array.reshape(shape) if shape is not None else array


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _reference_inputs(reference: dict) -> dict:
    spatial = tuple(reference["shape"])
    classes = reference["classes"]
    return {
        "image_values": _float(reference["images_bits"], (*spatial, 1)),
        "initial_memberships": _float(reference["initial_memberships_bits"], (*spatial, classes)),
        "class_types": reference["class_types"],
        "centroids": _float(reference["centroids_bits"], (1, classes)),
        "modality_weights": _float([reference["modality_weight_bits"]]),
        "class_variance_weights": _float(reference["variance_bits"], (1, classes)),
        "image_variances": _float([reference["image_variance_bits"]]),
        "transformed_shapes": _float(reference["shape_bits"], (*spatial, classes)),
        "closest_labels": np.asarray(reference["closest_labels"], dtype=np.int8).reshape((*spatial, 4)),
        "smoothing_priors": _float(reference["smoothing_priors_bits"], (classes,)),
        "centroid_similarity": _float(reference["similarity_bits"], (classes, classes)),
        "atlas_shape_priors": _float(reference["atlas_priors_bits"], (classes,)),
        "smoothing": reference["smoothing"],
        "prior_factor": reference["prior_factor"],
        "prior_coefficient": reference["prior_coefficient"],
        "initial_outlier_membership": _float(reference["initial_outlier_bits"], spatial),
        "outlier_value": reference["outlier_value"],
        "intensity_field": _float(reference["field_bits"], (*spatial, 1)),
        "psf_bias": _float(reference["bias_bits"], (*spatial, 1)),
    }


def test_supported_mixed_types_match_pinned_vandycruise_java_excerpt() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == SEGMENTATION_SHA256
    assert reference["class_types"] == ["OBJECT", "MASK", "OUTLIER"]
    inputs = _reference_inputs(reference)
    initial_before = inputs["initial_memberships"].copy()
    outlier_before = inputs["initial_outlier_membership"].copy()

    result = toads_exact_membership_dispatch(**inputs)

    assert _bits(result.memberships) == reference["updated_bits"]
    assert _bits(result.outlier_membership) == reference["updated_outlier_bits"]
    assert _bits(result.objective) == [reference["objective_bits"]]
    np.testing.assert_array_equal(inputs["initial_memberships"], initial_before)
    np.testing.assert_array_equal(inputs["initial_outlier_membership"], outlier_before)
    # Source iteration visits only x/y/z interior: both class and auxiliary
    # maps retain the input values on all six volume faces.
    np.testing.assert_array_equal(result.memberships[[0, -1]], initial_before[[0, -1]])
    np.testing.assert_array_equal(result.memberships[:, [0, -1]], initial_before[:, [0, -1]])
    np.testing.assert_array_equal(result.memberships[:, :, [0, -1]], initial_before[:, :, [0, -1]])
    np.testing.assert_array_equal(result.outlier_membership[[0, -1]], outlier_before[[0, -1]])
    np.testing.assert_array_equal(result.outlier_membership[:, [0, -1]], outlier_before[:, [0, -1]])
    np.testing.assert_array_equal(result.outlier_membership[:, :, [0, -1]], outlier_before[:, :, [0, -1]])


def _folded_inputs() -> tuple[dict, np.ndarray]:
    image_nii = nib.load(str(FOLDED / "images.nii.gz"))
    ribbon_nii = nib.load(str(FOLDED / "folded_ribbon.nii.gz"))
    assert nib.aff2axcodes(image_nii.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(ribbon_nii.affine) == ("R", "A", "S")
    images = np.asarray(image_nii.dataobj, dtype=np.float32)
    ribbon = np.asarray(ribbon_nii.dataobj, dtype=np.uint8)
    assert images.shape == (*ribbon.shape, 2)
    assert ribbon.sum() > 100
    # Five discrete lobes in a non-centrosymmetric shell exercise a spatially
    # varying fixed-input dispatch without invoking the full TOADS loop.
    assert np.unique(ribbon[7, :, :], return_counts=True)[0].size == 2
    assert not np.array_equal(ribbon, ribbon[::-1, :, :])

    nx, ny, nz = ribbon.shape
    classes = 2
    initial = np.empty((nx, ny, nz, classes), dtype=np.float32)
    shapes = np.empty_like(initial)
    labels = np.full((nx, ny, nz, 4), -1, dtype=np.int8)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                initial[x, y, z] = (0.2 + ((x + y + z) % 5) * 0.03, 0.3 + ((x + 2 * y + z) % 7) * 0.02)
                shell_value = np.float32(0.25 + 0.5 * ribbon[x, y, z])
                shapes[x, y, z] = (shell_value, np.float32(1.0) - shell_value)
                labels[x, y, z, 0] = (x + y + z) % 2
                if (x + 2 * y + z) % 3 == 0:
                    labels[x, y, z, 2] = 1 - labels[x, y, z, 0]

    return {
        "image_values": images,
        "initial_memberships": initial,
        "class_types": ["object", "mask"],
        "centroids": np.asarray([[0.42, 0.57], [0.61, 0.28]], dtype=np.float32),
        "modality_weights": np.asarray([0.7, 0.3], dtype=np.float32),
        "class_variance_weights": np.asarray([[0.8, 0.9], [1.1, 0.6]], dtype=np.float32),
        "image_variances": np.asarray([0.75, 0.25], dtype=np.float32),
        "transformed_shapes": shapes,
        "closest_labels": labels,
        "smoothing_priors": np.asarray([0.6, 0.85], dtype=np.float32),
        "centroid_similarity": np.asarray([[0, 0.4], [0.5, 0]], dtype=np.float32),
        "atlas_shape_priors": np.asarray([0.45, 0.55], dtype=np.float32),
        "smoothing": 0.03,
        "prior_factor": 0.62,
        "prior_coefficient": 0.08,
    }, initial


def test_object_mask_dispatch_composes_validated_kernel_on_folded_ras_fixture() -> None:
    inputs, initial = _folded_inputs()
    expected, expected_objective = toads_object_mask_membership_pass(
        inputs["image_values"],
        initial,
        inputs["centroids"],
        inputs["modality_weights"],
        inputs["class_variance_weights"],
        inputs["image_variances"],
        inputs["transformed_shapes"],
        inputs["closest_labels"],
        inputs["smoothing_priors"],
        inputs["centroid_similarity"],
        inputs["atlas_shape_priors"],
        inputs["smoothing"],
        inputs["prior_factor"],
        inputs["prior_coefficient"],
        class_types=("OBJECT", "MASK"),
    )
    result = toads_exact_membership_dispatch(**inputs)
    np.testing.assert_array_equal(result.memberships, expected)
    assert result.objective.view(np.uint32) == expected_objective.view(np.uint32)
    assert result.outlier_membership is None
    # Borders remain untouched for this ordinary two-type route as well.
    np.testing.assert_array_equal(result.memberships[0], initial[0])
    np.testing.assert_array_equal(result.memberships[-1], initial[-1])


@pytest.mark.parametrize(
    "kind",
    ["LESION", "OUTMASK", "RAYLEIGH", "OPTIMIZED_CSF", "OPTIMIZED_GM", "OPTIMIZED_WM"],
)
def test_uncomposed_source_branches_fail_explicitly(kind: str) -> None:
    inputs, _ = _folded_inputs()
    inputs["class_types"] = [kind, "OBJECT"]
    with pytest.raises(ValueError, match="unsupported exact-membership class types"):
        toads_exact_membership_dispatch(**inputs)


def test_outlier_auxiliary_arguments_must_match_dispatch_type() -> None:
    inputs, _ = _folded_inputs()
    inputs["initial_outlier_membership"] = np.zeros(inputs["image_values"].shape[:3], dtype=np.float32)
    inputs["outlier_value"] = 1.0
    with pytest.raises(ValueError, match="only valid when class_types contains OUTLIER"):
        toads_exact_membership_dispatch(**inputs)
