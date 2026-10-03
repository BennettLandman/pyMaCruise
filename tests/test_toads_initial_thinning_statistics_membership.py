from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_initial_thinning_composition import (
    toads_initial_fielded_membership_to_thinning,
)
from pymacruise.toads_initial_thinning_statistics_membership import (
    continue_initial_thinning_through_memberships,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_initial_thinning_startup"
JAVA_HARNESS = Path(__file__).parent / "java_reference" / "ToadsInitialStatisticsMembershipContinuationReferenceHarness.java"
JAVA_HARNESS_SHA256 = "89346fedd60534c8314d59ece0747a483649ef4b181efb408c7be52bbbf1360c"


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / name))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert image.header.get_xyzt_units()[0] == "mm"
    return np.asarray(image.dataobj)


def _sha(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f4").tobytes(order="C")).hexdigest()


def _java_float_digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype=">f4").tobytes(order="C")).hexdigest()


def _float_bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def test_startup_first_thinning_statistics_and_next_membership_compose() -> None:
    assert hashlib.sha256(JAVA_HARNESS.read_bytes()).hexdigest() == JAVA_HARNESS_SHA256
    spec = json.loads((ROOT / "spec.json").read_text())
    params = spec["parameters"]
    images = _load("images_ras.nii.gz").astype(np.float32)
    labels = _load("segmentation_ras.nii.gz").astype(np.uint8)
    initial = _load("initial_memberships_ras.nii.gz").astype(np.float32)
    fields = _load("fields_ras.nii.gz").astype(np.float32)
    bias = _load("psf_bias_ras.nii.gz").astype(np.float32)
    shapes = _load("transformed_shapes_ras.nii.gz").astype(np.float32)
    centers = np.asarray(params["centroids"], dtype=np.float32)
    exponents = np.asarray(params["exponent_by_channel_class"], dtype=np.float32)
    fixed = [images.copy(), labels.copy(), initial.copy(), fields.copy(), shapes.copy()]
    startup = toads_initial_fielded_membership_to_thinning(
        image_values=images,
        segmentation=labels,
        initial_memberships=initial,
        centroids=centers,
        exponent_by_channel_class=exponents,
        modality_weights=params["modality_weights"],
        class_variance_weights=params["class_variance_weights"],
        image_variances=params["image_variances"],
        transformed_shapes=shapes,
        smoothing_priors=params["smoothing_priors"],
        centroid_similarity=params["centroid_similarity"],
        atlas_shape_priors=params["atlas_shape_priors"],
        smoothing=spec["membership_smoothing"],
        prior_factor=spec["prior_factor"],
        prior_coefficient=spec["prior_coefficient"],
        channel_fields=[fields[..., 0], fields[..., 1]],
        limit=spec["relation_limit"],
        volume_ratio=spec["thinning_volume_ratio"],
        thinning_limit=spec["thinning_limit"],
        psf_bias=bias,
    )
    first_stage_reference = json.loads((ROOT / "java_reference.json").read_text())
    assert _java_float_digest(startup.membership_pass.memberships) == first_stage_reference["membership_sha256"]
    assert hashlib.sha256(startup.thinning.segmentation.astype(np.int8).tobytes(order="C")).hexdigest() == first_stage_reference["segmentation_sha256"]
    category_exponents = np.asarray(
        [exponents[:, 1], exponents[:, 3], exponents[:, 8]], dtype=np.float32
    ).T
    continued = continue_initial_thinning_through_memberships(
        startup,
        image_values=images,
        initial_centroids=centers,
        category_exponents=category_exponents,
        transformed_shapes=shapes,
        image_variance_weights=params["image_variances"],
        class_variance_weights=params["class_variance_weights"],
        smoothing_priors=params["smoothing_priors"],
        atlas_shape_priors=params["atlas_shape_priors"],
        exponent_by_channel_class=exponents,
        modality_weights=params["modality_weights"],
        prior_factor=spec["prior_factor"],
        prior_coefficient=spec["prior_coefficient"],
        smoothing=spec["membership_smoothing"],
        prior_scale=0.05,
        intensity_field=fields,
        psf_bias=bias,
    )
    java_reference = json.loads((ROOT / "continuation_java_reference.json").read_text())

    assert continued.stage_order == (
        *startup.stage_order,
        "computeExactCentroids",
        "computeExactVariances",
        "computeCentroidSimilarity",
        "computeExactMemberships",
    )
    assert continued.centroids.shape == continued.class_standard_deviations.shape == (2, 11)
    assert continued.centroid_similarity.shape == (11, 11)
    assert continued.membership_pass.memberships.shape == initial.shape
    assert np.isfinite(continued.centroids).all()
    assert np.isfinite(continued.class_standard_deviations).all()
    assert np.isfinite(continued.membership_pass.memberships).all()
    assert java_reference["source_order"] == [
        "P18cyn:startup membership and first thinning",
        "computeExactCentroids",
        "computeExactVariances",
        "computeCentroidSimilarity",
        "computeExactMemberships",
    ]
    assert java_reference["source_sha256"] == first_stage_reference["toads_source_sha256"]
    assert java_reference["startup_membership_sha256"] == first_stage_reference["membership_sha256"]
    assert java_reference["post_thinning_segmentation_sha256"] == first_stage_reference["segmentation_sha256"]
    assert _float_bits(continued.centroids) == [
        value for channel in java_reference["updated_centroids_bits"] for value in channel
    ]
    assert _float_bits(continued.class_standard_deviations) == [
        value for channel in java_reference["standard_deviation_bits"] for value in channel
    ]
    assert _java_float_digest(continued.centroid_similarity) == java_reference["centroid_similarity_sha256"]
    assert _java_float_digest(continued.membership_pass.class_data_energy) == java_reference["class_energy_sha256"]
    assert _java_float_digest(continued.membership_pass.denominator) == java_reference["denominator_sha256"]
    assert _java_float_digest(continued.membership_pass.memberships) == java_reference["updated_membership_sha256"]
    assert int(np.asarray(continued.membership_pass.objective, dtype=np.float32).view(np.int32)) == java_reference["objective_bits"]
    np.testing.assert_allclose(
        continued.membership_pass.memberships[1:-1, 1:-1, 1:-1].sum(axis=-1),
        1.0,
        rtol=0,
        atol=3e-7,
    )
    for current, before in zip((images, labels, initial, fields, shapes), fixed, strict=True):
        np.testing.assert_array_equal(current, before)
    assert len(set(_sha(value) for value in (continued.centroids, continued.class_standard_deviations))) == 2


def test_continuation_requires_the_selected_startup_result() -> None:
    with np.testing.assert_raises_regex(TypeError, "selected startup composition"):
        continue_initial_thinning_through_memberships(
            object(),
            image_values=np.zeros((3, 3, 3, 2), dtype=np.float32),
            initial_centroids=np.zeros((2, 11), dtype=np.float32),
            category_exponents=np.ones((2, 3), dtype=np.float32),
            transformed_shapes=np.zeros((3, 3, 3, 11), dtype=np.float32),
            image_variance_weights=np.ones(2, dtype=np.float32),
            class_variance_weights=np.ones((2, 11), dtype=np.float32),
            smoothing_priors=np.ones(11, dtype=np.float32),
            atlas_shape_priors=np.ones(11, dtype=np.float32),
            exponent_by_channel_class=np.ones((2, 11), dtype=np.float32),
            modality_weights=np.ones(2, dtype=np.float32),
            prior_factor=1.0,
            prior_coefficient=0.0,
            smoothing=0.0,
            prior_scale=0.05,
        )
