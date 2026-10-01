from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_approximate_centroids import update_default_toads_approximate_centroids
from pymacruise.toads_approximate_variances import compute_default_toads_approximate_variances


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_approximate_variances"
CENTROID_ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_approximate_centroids"
JAVA = json.loads((ROOT / "java_reference.json").read_text())
CENTROID_JAVA = json.loads((CENTROID_ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "ToadsDefaultApproximateVarianceReferenceHarness.java"
HARNESS_SHA256 = "aad6815d6812534b41223a0430804799c2c14ff6b840168f1e23256644d9eefd"
SOURCE_HASHES = {
    "segmentation_sha256": "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884",
    "algorithm_sha256": "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d",
    "wrapper_sha256": "d00a7ad609ab0c6381334248b53f4daadc070161af67ef853e2a552dd008d115",
    "descriptor_sha256": "cfe3e2556cca554122f5f16427878b6924cb1de6947e9392b2918bed24a4a939",
}


def _load(name: str, *, dtype=None) -> np.ndarray:
    volume = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(volume.affine) == ("R", "A", "S")
    np.testing.assert_allclose(
        nib.affines.voxel_sizes(volume.affine), (.85, 1.15, 1.55), rtol=0, atol=1e-6
    )
    array = np.asarray(volume.dataobj)
    return array.astype(dtype, copy=False) if dtype is not None else array.copy()


def _bits(values: list[int]) -> np.ndarray:
    return np.asarray(values, dtype=np.int32).view(np.float32).reshape(2, 11)


def test_selected_approximate_variances_match_java_reference_bits() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    for key, expected in SOURCE_HASHES.items():
        assert JAVA[key] == expected
    assert JAVA["shape"] == [18, 18, 18, 2]
    assert JAVA["classes"] == 11 and JAVA["approx"] == 4
    assert JAVA["epsilon"] == 1e-30

    images = _load("images", dtype=np.float32)
    field = _load("field", dtype=np.float32)
    memberships = _load("memberships", dtype=np.float32)
    best = _load("best")
    segmentation = _load("segmentation")
    centroids = np.loadtxt(ROOT / "updated_centroids.txt", dtype=np.float32)
    old_centroids = np.loadtxt(CENTROID_ROOT / "centroids.txt", dtype=np.float32)
    exponents = np.loadtxt(CENTROID_ROOT / "optimized_exponents.txt", dtype=np.float32)
    centroid_from_python = update_default_toads_approximate_centroids(
        images, memberships, best, old_centroids, field, exponents
    ).centroids
    expected_centroids = _bits(CENTROID_JAVA["centroid_bits"])
    np.testing.assert_array_equal(centroid_from_python.view(np.int32), expected_centroids.view(np.int32))
    np.testing.assert_array_equal(centroids.view(np.int32), expected_centroids.view(np.int32))
    originals = tuple(x.copy() for x in (images, field, memberships, best, segmentation, centroids))

    result = compute_default_toads_approximate_variances(
        images, memberships, best, segmentation, centroids, intensity_field=field
    )
    np.testing.assert_array_equal(result.numerators.view(np.int32), _bits(JAVA["field_numerator_bits"]).view(np.int32))
    np.testing.assert_array_equal(result.denominators.view(np.int32), _bits(JAVA["field_denominator_bits"]).view(np.int32))
    np.testing.assert_array_equal(result.standard_deviations.view(np.int32), _bits(JAVA["field_stddev_bits"]).view(np.int32))

    raw_result = compute_default_toads_approximate_variances(
        images, memberships, best, segmentation, centroids
    )
    np.testing.assert_array_equal(raw_result.standard_deviations.view(np.int32), _bits(JAVA["raw_stddev_bits"]).view(np.int32))
    assert not np.array_equal(raw_result.standard_deviations, result.standard_deviations)

    matching = best == segmentation[..., None]
    expected_matches = int(matching.sum()) * images.shape[3]
    boundary = np.zeros(segmentation.shape, dtype=bool)
    boundary[[0, -1], :, :] = True
    boundary[:, [0, -1], :] = True
    boundary[:, :, [0, -1]] = True
    boundary_matches = int(np.count_nonzero(matching & boundary[..., None])) * images.shape[3]
    duplicate_matches = int(np.count_nonzero(np.sum(matching, axis=-1) > 1)) * images.shape[3]
    assert JAVA["matching_slots"] == expected_matches
    assert JAVA["mismatching_slots"] == int(best.size * images.shape[3]) - expected_matches
    assert JAVA["boundary_matching_slots"] == boundary_matches > 0
    assert JAVA["duplicate_matching_voxels"] == duplicate_matches > 0

    # In the Java edge pass, every matching MASK slot has zero membership.
    zero_memberships = memberships.copy()
    zero_memberships[(best == 0) & (segmentation[..., None] == 0)] = np.float32(0.0)
    zero_result = compute_default_toads_approximate_variances(
        images, zero_memberships, best, segmentation, centroids, intensity_field=field
    )
    np.testing.assert_array_equal(
        zero_result.denominators.view(np.int32), _bits(JAVA["zero_membership_denominator_bits"]).view(np.int32)
    )
    np.testing.assert_array_equal(
        zero_result.standard_deviations.view(np.int32), _bits(JAVA["zero_membership_stddev_bits"]).view(np.int32)
    )
    assert np.all(zero_result.denominators[:, 0] > np.float32(0.0))

    # Relabel all class-10 voxels to class 9 so no candidate matches class 10.
    no_support_labels = segmentation.copy()
    no_support_labels[no_support_labels == 10] = 9
    no_support = compute_default_toads_approximate_variances(
        images, memberships, best, no_support_labels, centroids, intensity_field=field
    )
    assert no_support.standard_deviations[0, 10].view(np.uint32) == JAVA["no_support_class10_stddev_bits"][0]
    assert no_support.standard_deviations[1, 10].view(np.uint32) == JAVA["no_support_class10_stddev_bits"][1]
    assert np.array_equal(no_support.standard_deviations[:, 10], np.zeros(2, dtype=np.float32))

    for current, original in zip(
        (images, field, memberships, best, segmentation, centroids), originals, strict=True
    ):
        np.testing.assert_array_equal(current, original)
