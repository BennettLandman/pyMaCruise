from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_default_atlas_membership import (
    DEFAULT_ATLAS_CLASS_TYPES,
    toads_default_atlas_membership_pass,
)
from pymacruise.toads_default_exact_centroids import toads_default_exact_centroid_update
from pymacruise.toads_default_exact_variances import toads_default_exact_class_variances
from pymacruise.toads_math import toads_optimized_centroid_similarity


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_11_continue_thinning"
REFERENCE = json.loads((ROOT / "membership_reference.json").read_text())
CHAIN_REFERENCE = json.loads((ROOT / "java_reference.json").read_text())


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / name))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    return np.asarray(image.dataobj)


def _sha(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f4").tobytes(order="C")).hexdigest()


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _labels(values: str) -> np.ndarray:
    decoded = np.empty(len(values), dtype=np.int16)
    for i, value in enumerate(values):
        decoded[i] = ord(value) - ord("0") if value <= "9" else ord(value) - ord("a") + 10
    return decoded.reshape((14, 14, 14))


def test_p18cym_two_channel_membership_refresh_matches_java_reference() -> None:
    assert REFERENCE["source_order"] == [
        "P18cyl:propagateThinning",
        "computeExactCentroids",
        "computeExactVariances",
        "fixed:second synthetic channel and exact pass parameters",
        "computeExactMemberships",
    ]
    assert REFERENCE["source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert REFERENCE["power_table_source_sha256"] == "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
    harness = Path(__file__).parent / "java_reference" / "ToadsDefault11MembershipContinuationReferenceHarness.java"
    assert REFERENCE["harness_sha256"] == hashlib.sha256(harness.read_bytes()).hexdigest()
    assert REFERENCE["post_thinning_segmentation"] == CHAIN_REFERENCE["segmentation"]
    assert REFERENCE["updated_centroids_channel0_bits"] == CHAIN_REFERENCE["updated_centroids_bits"]
    assert REFERENCE["standard_deviation_channel0_bits"] == CHAIN_REFERENCE["standard_deviation_bits"]

    images = np.stack(
        (_load("image_values.nii.gz"), _load("image_values_channel2.nii.gz")), axis=-1
    ).astype(np.float32)
    field = np.stack(
        (_load("intensity_field.nii.gz"), _load("intensity_field_channel2.nii.gz")), axis=-1
    ).astype(np.float32)
    initial_memberships = _load("exact_memberships.nii.gz").astype(np.float32)
    shapes = _load("transformed_shapes.nii.gz").astype(np.float32)
    closest = _load("membership_closest_labels.nii.gz").astype(np.int8)
    segmentation = _labels(CHAIN_REFERENCE["segmentation"])
    assert images.shape == field.shape == (14, 14, 14, 2)
    assert initial_memberships.shape == shapes.shape == (14, 14, 14, 11)
    assert len(DEFAULT_ATLAS_CLASS_TYPES) == 11
    assert closest.shape == (14, 14, 14, 4)
    assert np.all(shapes == 0)
    assert np.allclose(initial_memberships.sum(axis=-1), 1.0, rtol=0, atol=2e-7)
    assert _sha(images[..., 0]) == REFERENCE["image_channel1_sha256"]
    assert _sha(field[..., 0]) == REFERENCE["field_channel1_sha256"]
    assert _sha(images) == REFERENCE["image_channel2_sha256"]
    assert _sha(field) == REFERENCE["field_channel2_sha256"]
    assert _sha(initial_memberships) == REFERENCE["initial_membership_sha256"]
    assert _sha(shapes) == REFERENCE["shape_sha256"]

    params = json.loads((ROOT / "membership_parameters.json").read_text())
    stats = json.loads((ROOT / "statistics_inputs.json").read_text())
    centers_before = np.stack(
        (
            np.asarray(stats["initial_centroids"], dtype=np.float32),
            np.asarray([(11 - k) / 12.0 for k in range(11)], dtype=np.float32),
        )
    )
    category = np.asarray(params["category_exponents"], dtype=np.float32)
    exponents = np.empty((2, 11), dtype=np.float32)
    exponents[:, 0] = 1.0
    exponents[:, 1:3] = category[0]
    exponents[:, 3:8] = category[1]
    exponents[:, 8:11] = category[2]
    source_inputs = [images.copy(), field.copy(), initial_memberships.copy(), shapes.copy(), closest.copy(), segmentation.copy()]

    centers = toads_default_exact_centroid_update(
        images, initial_memberships, centers_before, category[None, :].repeat(2, axis=0),
        intensity_field=field, smoothness=stats["centroid_smoothness"],
    )
    deviations = toads_default_exact_class_variances(
        images, initial_memberships, centers, segmentation, intensity_field=field,
    )
    assert _bits(centers[0]) == REFERENCE["updated_centroids_channel0_bits"]
    assert _bits(centers[1]) == REFERENCE["updated_centroids_channel1_bits"]
    assert _bits(deviations[0]) == REFERENCE["standard_deviation_channel0_bits"]
    assert _bits(deviations[1]) == REFERENCE["standard_deviation_channel1_bits"]

    similarity = toads_optimized_centroid_similarity(
        centers, deviations, np.asarray(params["image_variances"], dtype=np.float32), params["prior_scale"]
    )
    result = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=initial_memberships,
        centroids=centers,
        exponent_by_channel_class=exponents,
        modality_weights=np.asarray(params["modality_weights"], dtype=np.float32),
        class_variance_weights=np.asarray(params["class_variance_weights"], dtype=np.float32),
        image_variances=np.asarray(params["image_variances"], dtype=np.float32),
        transformed_shapes=shapes,
        closest_labels=closest,
        smoothing_priors=np.asarray(params["smoothing_priors"], dtype=np.float32),
        centroid_similarity=similarity,
        atlas_shape_priors=np.asarray(params["atlas_shape_priors"], dtype=np.float32),
        smoothing=params["smoothing"],
        prior_factor=params["prior_factor"],
        prior_coefficient=params["prior_coefficient"],
        intensity_field=field,
    )
    assert _sha(result.memberships) == REFERENCE["membership_bits_sha256"]
    assert _sha(result.denominator) == REFERENCE["denominator_sha256"]
    assert _sha(result.class_data_energy) == REFERENCE["class_energy_sha256"]
    assert int(np.asarray(result.objective, dtype=np.float32).view(np.int32)) == REFERENCE["objective_bits"]
    assert result.memberships.dtype == result.class_data_energy.dtype == result.denominator.dtype == np.float32
    assert not np.array_equal(result.memberships[1:-1, 1:-1, 1:-1], initial_memberships[1:-1, 1:-1, 1:-1])
    np.testing.assert_array_equal(result.memberships[0], initial_memberships[0])
    np.testing.assert_array_equal(result.memberships[-1], initial_memberships[-1])
    for current, fixed in zip((images, field, initial_memberships, shapes, closest, segmentation), source_inputs, strict=True):
        np.testing.assert_array_equal(current, fixed)
