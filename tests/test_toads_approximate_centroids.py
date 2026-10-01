from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_approximate_centroids import (
    DEFAULT_TOADS_11_CLASS_TYPES,
    update_default_toads_approximate_centroids,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_approximate_centroids"
JAVA = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "ToadsDefaultApproximateCentroidReferenceHarness.java"
HARNESS_SHA256 = "a63cd48f70ccc369821513878abeb421a1486a1bb738c0d11a5d4f2afc85cb85"
P18CYB = Path(__file__).parent / "data" / "synthetic" / "toads_inhomogeneity_transition" / "java_reference.json"
EXPECTED_SOURCE_HASHES = {
    "segmentation_sha256": "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884",
    "algorithm_sha256": "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d",
    "wrapper_sha256": "d00a7ad609ab0c6381334248b53f4daadc070161af67ef853e2a552dd008d115",
    "power_vandy_sha256": "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678",
    "power_copy_sha256": "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae",
    "descriptor_sha256": "cfe3e2556cca554122f5f16427878b6924cb1de6947e9392b2918bed24a4a939",
}


def _load(name: str, *, dtype=None) -> np.ndarray:
    volume = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(volume.affine) == ("R", "A", "S")
    np.testing.assert_allclose(
        nib.affines.voxel_sizes(volume.affine), (.85, 1.15, 1.55), rtol=0, atol=1e-6
    )
    values = np.asarray(volume.dataobj)
    return values.astype(dtype, copy=False) if dtype is not None else values.copy()


def _float_bits(values: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(values, dtype=np.int32).view(np.float32).reshape(shape)


def test_selected_approximate_centroids_match_java_reference() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    for key, expected in EXPECTED_SOURCE_HASHES.items():
        assert JAVA[key] == expected
    assert JAVA["shape"] == [18, 18, 18, 2]
    assert JAVA["classes"] == 11 and JAVA["approx"] == 4
    assert JAVA["centroid_mode"] == "damped" and JAVA["smoothness"] == 0.1
    assert DEFAULT_TOADS_11_CLASS_TYPES == (
        "MASK", "OPTIMIZED_CSF", "OPTIMIZED_CSF", "OPTIMIZED_GM", "OPTIMIZED_GM",
        "OPTIMIZED_GM", "OPTIMIZED_GM", "OPTIMIZED_GM", "OPTIMIZED_WM",
        "OPTIMIZED_WM", "OPTIMIZED_WM",
    )

    images = _load("images", dtype=np.float32)
    field = _load("field", dtype=np.float32)
    memberships = _load("memberships", dtype=np.float32)
    best = _load("best")
    segmentation = _load("segmentation")
    centroids = np.loadtxt(ROOT / "centroids.txt", dtype=np.float32)
    exponents = np.loadtxt(ROOT / "optimized_exponents.txt", dtype=np.float32)
    source = json.loads(P18CYB.read_text())
    for channel in range(2):
        digest = hashlib.sha256(np.asarray(field[..., channel], dtype="<f4").tobytes(order="F")).hexdigest()
        assert digest == source["final_field_sha256"][channel]
    assert np.any(field != np.float32(1.0))
    assert not np.array_equal(field[..., 0], field[..., 1])
    originals = tuple(x.copy() for x in (images, memberships, best, segmentation, field, centroids))

    result = update_default_toads_approximate_centroids(
        images, memberships, best, centroids, field, exponents,
    )
    expected_shape = (2, 11)
    np.testing.assert_array_equal(
        result.numerators.view(np.int32), _float_bits(JAVA["numerator_bits"], expected_shape).view(np.int32)
    )
    np.testing.assert_array_equal(
        result.denominators.view(np.int32), _float_bits(JAVA["denominator_bits"], expected_shape).view(np.int32)
    )
    np.testing.assert_array_equal(
        result.centroids.view(np.int32), _float_bits(JAVA["centroid_bits"], expected_shape).view(np.int32)
    )
    assert JAVA["zero_fallbacks"] == 2
    assert JAVA["table_lookups"] > 0
    assert JAVA["boundary_slots"] > 0
    assert JAVA["duplicate_candidate_voxels"] > 0

    # The Java reference explicitly suppresses class 10's candidate weights;
    # with no denominator its old centroid must remain exactly unchanged.
    no_class10_support = memberships.copy()
    no_class10_support[best == 10] = np.float32(0.0)
    empty = update_default_toads_approximate_centroids(
        images, no_class10_support, best, centroids, field, exponents,
    )
    assert empty.denominators[0, 10] == np.float32(0.0)
    assert empty.centroids[0, 10].view(np.uint32) == JAVA["empty_class10_centroid_bits"]
    assert empty.centroids[1, 10].view(np.uint32) == JAVA["empty_class10_centroid_bits_ch1"]

    # The folded fixture has exactly two optimized-category distance samples
    # at or below source ZERO; the remaining candidates exercise PowerTable.
    optimized_classes = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
    zero_count = 0
    for c in range(images.shape[3]):
        for x in range(images.shape[0]):
            for y in range(images.shape[1]):
                for z in range(images.shape[2]):
                    corrected = np.float32(field[x, y, z, c] * images[x, y, z, c])
                    for slot in range(4):
                        klass = int(best[x, y, z, slot])
                        if klass in optimized_classes:
                            residual = np.float32(corrected - centroids[c, klass])
                            distance = np.float32(residual * residual)
                            zero_count += int(distance <= np.float32(1e-10))
    assert zero_count == JAVA["zero_fallbacks"]
    for current, original in zip(
        (images, memberships, best, segmentation, field, centroids), originals, strict=True
    ):
        np.testing.assert_array_equal(current, original)


def test_approximate_centroid_rejects_bad_candidate_shapes_and_classes() -> None:
    images = _load("images", dtype=np.float32)
    field = _load("field", dtype=np.float32)
    memberships = _load("memberships", dtype=np.float32)
    best = _load("best")
    centroids = np.loadtxt(ROOT / "centroids.txt", dtype=np.float32)
    exponents = np.loadtxt(ROOT / "optimized_exponents.txt", dtype=np.float32)
    import pytest

    with pytest.raises(ValueError, match="four candidate slots"):
        update_default_toads_approximate_centroids(
            images, memberships[..., :3], best[..., :3], centroids, field, exponents
        )
    invalid = best.copy()
    invalid[0, 0, 0, 0] = 11
    with pytest.raises(ValueError, match="outside the selected"):
        update_default_toads_approximate_centroids(
            images, memberships, invalid, centroids, field, exponents
        )
