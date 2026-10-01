from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_inhomogeneity_fit_degree23 import (
    fit_toads_approximate_chebyshev_field_degree23,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_inhomogeneity_fit_degree23"
JAVA = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "ToadsInhomogeneityDegree23FitReferenceHarness.java"
HARNESS_SHA256 = "2b88fbdd93d69eb5ad8fd44bd80c9abd5cc5ff8eb43157c9c99ffff111fee9d9"
COEFFICIENT_MAX_ABS = {2: 3.0e-13, 3: 4.0e-11}
FIELD_MAX_ULP = {2: 0, 3: 0}
SOURCE_HASHES = {
    "source_sha256": "6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8",
    "jama_wrapper_sha256": "efd52d4a8f86984c1dc52d13e01175e07120ccc1e9867283846a0b17f1025526",
    "jama_matrix_sha256": "76ab1241a91518c3b71db14d19a2a3d6c35376f65fdc9804fe26fa79e964ff9c",
    "jama_lu_sha256": "9b41fefaef78fdc95a2c8fdfebb1aaee8375dcaf6c7b8eddb30acdf7f5a077a6",
    "power_table_sha256": "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678",
}


def _load_inputs():
    names = ("image", "memberships", "best", "segmentation")
    volumes = tuple(nib.load(str(ROOT / f"{name}.nii.gz")) for name in names)
    for volume in volumes:
        assert nib.aff2axcodes(volume.affine) == ("R", "A", "S")
        np.testing.assert_allclose(
            nib.affines.voxel_sizes(volume.affine), (.85, 1.15, 1.55), rtol=0, atol=1e-6
        )
    return tuple(np.asarray(volume.dataobj) for volume in volumes)


def _from_double_bits(values: list[int]) -> np.ndarray:
    return np.asarray(values, dtype=np.int64).view(np.float64)


def _from_java_float_bits(values: list[int], shape: tuple[int, int, int]) -> np.ndarray:
    return np.asarray(values, dtype=np.int32).view(np.float32).reshape(shape, order="F")


@pytest.mark.parametrize("degree", [2, 3])
def test_degree23_field_fit_matches_selected_java_reference(degree: int) -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    for key, expected in SOURCE_HASHES.items():
        assert JAVA[key] == expected
    assert JAVA["shape"] == [18, 18, 18]
    image, memberships, best, segmentation = _load_inputs()
    centroids = np.loadtxt(ROOT / "centroids.txt", dtype=np.float32)
    original = tuple(value.copy() for value in (image, memberships, best, segmentation, centroids))

    actual = fit_toads_approximate_chebyshev_field_degree23(
        image, memberships, best, segmentation, centroids, degree=degree,
    )
    reference = next(entry for entry in JAVA["degrees"] if entry["degree"] == degree)
    n_parameters = 10 if degree == 2 else 20
    assert reference["parameter_count"] == n_parameters
    assert int(actual.correction_mask.sum()) == reference["mask_count"]
    assert actual.right_hand_side.shape == (n_parameters,)
    assert actual.normal_matrix.shape == (n_parameters, n_parameters)
    np.testing.assert_array_equal(actual.right_hand_side, _from_double_bits(reference["rhs_bits"]))
    np.testing.assert_array_equal(
        actual.normal_matrix.ravel(), _from_double_bits(reference["normal_matrix_bits"])
    )

    expected_coefficients = _from_double_bits(reference["coefficient_bits"])
    coefficient_max_abs = float(np.max(np.abs(actual.coefficients - expected_coefficients)))
    assert coefficient_max_abs <= COEFFICIENT_MAX_ABS[degree]

    expected_field = _from_java_float_bits(reference["field_bits"], image.shape)
    actual_bits = actual.field.view(np.uint32).astype(np.int64)
    expected_bits = expected_field.view(np.uint32).astype(np.int64)
    field_ulp_max = int(np.max(np.abs(actual_bits - expected_bits)))
    assert field_ulp_max <= FIELD_MAX_ULP[degree]
    assert reference["sample_count"] == 94
    np.testing.assert_array_equal(actual.field[~actual.correction_mask], 1.0)
    assert abs(float(np.mean(actual.field[actual.correction_mask], dtype=np.float64)) - 1.0) < 2e-6
    for actual_input, original_input in zip(
        (image, memberships, best, segmentation, centroids), original, strict=True
    ):
        np.testing.assert_array_equal(actual_input, original_input)


def test_degree23_fit_rejects_singular_support_and_bad_class_indices() -> None:
    image, memberships, best, segmentation = _load_inputs()
    centroids = np.loadtxt(ROOT / "centroids.txt", dtype=np.float32)
    one_voxel = np.zeros_like(segmentation)
    one_voxel[0, 0, 0] = 1
    with pytest.raises(ValueError, match="singular"):
        fit_toads_approximate_chebyshev_field_degree23(
            image, memberships, best, one_voxel, centroids, degree=2,
        )
    invalid_best = best.copy()
    invalid_best[0, 0, 0, 1] = centroids.size
    with pytest.raises(ValueError, match="outside the centroid"):
        fit_toads_approximate_chebyshev_field_degree23(
            image, memberships, invalid_best, segmentation, centroids, degree=3,
        )


def test_degree23_fit_rejects_unsupported_degree_and_shape_mismatch() -> None:
    image, memberships, best, segmentation = _load_inputs()
    centroids = np.loadtxt(ROOT / "centroids.txt", dtype=np.float32)
    with pytest.raises(ValueError, match="degree must be 2 or 3"):
        fit_toads_approximate_chebyshev_field_degree23(
            image, memberships, best, segmentation, centroids, degree=1,
        )
    with pytest.raises(ValueError, match="match the image grid"):
        fit_toads_approximate_chebyshev_field_degree23(
            image, memberships[:-1], best[:-1], segmentation, centroids, degree=2,
        )
