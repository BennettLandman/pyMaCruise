"""Selected Java comparison for TOADS' degree-one approximate field fit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_inhomogeneity_fit import fit_toads_approximate_chebyshev_field_degree1


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_inhomogeneity_fit_degree1"
JAVA = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "ToadsInhomogeneityDegreeOneFitReferenceHarness.java"
HARNESS_SHA256 = "80b40a15ee803fba404aba312307b15a6b7c7d93c230d2c6432d7eab770c6910"
SOURCE_HASHES = {
    "source_sha256": "6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8",
    "jama_wrapper_sha256": "efd52d4a8f86984c1dc52d13e01175e07120ccc1e9867283846a0b17f1025526",
    "jama_matrix_sha256": "76ab1241a91518c3b71db14d19a2a3d6c35376f65fdc9804fe26fa79e964ff9c",
    "jama_lu_sha256": "9b41fefaef78fdc95a2c8fdfebb1aaee8375dcaf6c7b8eddb30acdf7f5a077a6",
    "power_table_sha256": "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678",
}


def _inputs():
    names = ("image", "memberships", "best", "segmentation")
    volumes = tuple(nib.load(str(ROOT / f"{name}.nii.gz")) for name in names)
    for volume in volumes:
        assert nib.aff2axcodes(volume.affine) == ("R", "A", "S")
        np.testing.assert_allclose(nib.affines.voxel_sizes(volume.affine), (.85, 1.15, 1.55), rtol=0, atol=1e-6)
    return tuple(np.asarray(volume.dataobj) for volume in volumes)


def _field_from_java_bits(values: list[int], shape: tuple[int, int, int]) -> np.ndarray:
    # Java's flat volume is indexed x + nx * (y + ny*z); restore x/y/z order.
    return np.asarray(values, dtype=np.int32).view(np.float32).reshape(shape, order="F")


def test_degree_one_field_fit_matches_selected_java_reference() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    for key, expected in SOURCE_HASHES.items():
        assert JAVA[key] == expected
    images, memberships, best, segmentation = _inputs()
    centroids = np.loadtxt(ROOT / "centroids.txt", dtype=np.float32)
    originals = tuple(value.copy() for value in (images, memberships, best, segmentation, centroids))
    result = fit_toads_approximate_chebyshev_field_degree1(
        images, memberships, best, segmentation, centroids,
    )

    assert JAVA["degree"] == 1
    assert JAVA["stride"] == 3
    assert int(result.correction_mask.sum()) == JAVA["mask_count"]
    expected_rhs = np.asarray(JAVA["rhs_bits"], dtype=np.int64).view(np.float64)
    expected_matrix = np.asarray(JAVA["normal_matrix_bits"], dtype=np.int64).view(np.float64)
    expected_coefficients = np.asarray(JAVA["coefficient_bits"], dtype=np.int64).view(np.float64)
    np.testing.assert_array_equal(result.right_hand_side, expected_rhs)
    np.testing.assert_array_equal(result.normal_matrix.ravel(), expected_matrix)
    np.testing.assert_allclose(result.coefficients, expected_coefficients, rtol=0, atol=8e-16)

    java_field = _field_from_java_bits(JAVA["field_bits"], images.shape)
    np.testing.assert_allclose(result.field, java_field, rtol=0, atol=1.2e-7)
    np.testing.assert_array_equal(result.field[~result.correction_mask], 1.0)
    assert abs(float(np.mean(result.field[result.correction_mask], dtype=np.float64)) - 1.0) < 7e-7
    # The all-grid field pass includes valid, non-subsampled points in its mean.
    xyz = np.indices(images.shape)
    unsampled_masked = result.correction_mask & np.any(xyz % 3 != 0, axis=0)
    assert np.any(unsampled_masked)
    assert np.any(result.field[unsampled_masked] != 1.0)
    # Slot zero is the atlas MASK class and contributes no fit weight.
    assert np.all(best[..., 0] == 0)
    assert np.all(best[..., 1:] > 0)
    for actual, original in zip((images, memberships, best, segmentation, centroids), originals, strict=True):
        np.testing.assert_array_equal(actual, original)


def test_degree_one_fit_rejects_empty_or_singular_background_support() -> None:
    image, memberships, best, segmentation = _inputs()
    centroids = np.loadtxt(ROOT / "centroids.txt", dtype=np.float32)
    with pytest.raises(ValueError, match="nonempty correction support"):
        fit_toads_approximate_chebyshev_field_degree1(
            image, memberships, best, np.zeros_like(segmentation), centroids,
        )
    no_samples = np.zeros_like(segmentation)
    no_samples[1, 1, 1] = 1
    with pytest.raises(ValueError, match="singular"):
        fit_toads_approximate_chebyshev_field_degree1(
            image, memberships, best, no_samples, centroids,
        )
    only_mask_candidates = np.zeros_like(best)
    with pytest.raises(ValueError, match="singular"):
        fit_toads_approximate_chebyshev_field_degree1(
            image, memberships, only_mask_candidates, segmentation, centroids,
        )


@pytest.mark.parametrize(
    ("image", "memberships", "best", "segmentation", "centroids", "message"),
    [
        (np.ones((3, 3)), np.ones((3, 3, 3, 2)), np.zeros((3, 3, 3, 2), np.uint8), np.ones((3, 3, 3), np.uint8), np.ones(2), "3D grid"),
        (np.ones((9, 9, 9)), np.ones((8, 9, 9, 2)), np.zeros((8, 9, 9, 2), np.uint8), np.ones((9, 9, 9), np.uint8), np.ones(3), "match the image grid"),
        (np.ones((9, 9, 9)), np.ones((9, 9, 9, 2)), np.full((9, 9, 9, 2), 4, np.uint8), np.ones((9, 9, 9), np.uint8), np.ones(3), "outside the centroid"),
    ],
)
def test_degree_one_fit_validates_inputs(image, memberships, best, segmentation, centroids, message) -> None:
    with pytest.raises(ValueError, match=message):
        fit_toads_approximate_chebyshev_field_degree1(image, memberships, best, segmentation, centroids)
