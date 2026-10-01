from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_default_atlas_membership import toads_default_atlas_membership_pass
from pymacruise.toads_default_exact_centroids import toads_default_exact_centroid_update
from pymacruise.toads_default_exact_variances import toads_default_exact_class_variances
from pymacruise.toads_math import toads_optimized_centroid_similarity, toads_optimized_image_variance
from pymacruise.toads_multichannel_exact_field_refit import toads_multichannel_exact_field_refit


ROOT = Path(__file__).parent
CASE_DIR = ROOT / "data/synthetic/toads_p18czo_nonzero_shapes"
CASE = json.loads((CASE_DIR / "case.json").read_text())
BASE_CASE_DIR = ROOT / "data/synthetic/toads_p18czo_multichannel_exact_field"
BASE_CASE = json.loads((BASE_CASE_DIR / "case.json").read_text())
JAVA_BITS = json.loads((BASE_CASE_DIR / "image_variance_reference.json").read_text())
SOURCE_DIR = ROOT / "data/synthetic" / BASE_CASE["source_fixture_directory"]
JAVA_HARNESS = ROOT.parent / "scripts/run_java_toads_p18czo_nonzero_shapes_multichannel_reference.sh"
SOURCE_SHA256 = "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
WRAPPER_SOURCE_SHA256 = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
POWER_TABLE_SHA256 = "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"


def _load(path: Path) -> nib.Nifti1Image:
    return nib.load(str(path))


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).astype(int).tolist()


def _sha(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f4").tobytes(order="C")).hexdigest()


def test_two_channel_dynamic_membership_uses_nonzero_ras_shapes_and_smoothing_priors(tmp_path: Path) -> None:
    channel_names = [entry.split(" / ")[0] for entry in BASE_CASE["channel_order"]]
    field_names = [entry.split(" / ")[1] for entry in BASE_CASE["channel_order"]]
    image_volumes = tuple(_load(SOURCE_DIR / name) for name in channel_names)
    field_volumes = tuple(_load(SOURCE_DIR / name) for name in field_names)
    membership_volume = _load(SOURCE_DIR / "exact_memberships.nii.gz")
    segmentation_volume = _load(BASE_CASE_DIR / BASE_CASE["segmentation"])
    shape_volume = _load(CASE_DIR / CASE["shape_file"])
    affine = image_volumes[0].affine
    for volume in (*image_volumes, *field_volumes, membership_volume, segmentation_volume, shape_volume):
        assert nib.aff2axcodes(volume.affine) == tuple(CASE["axis_codes"])
        np.testing.assert_allclose(volume.affine, affine, rtol=0, atol=1e-6)
    np.testing.assert_allclose(nib.affines.voxel_sizes(affine), CASE["voxel_spacing_mm"], atol=1e-6)

    images = np.stack(tuple(np.asarray(v.dataobj, dtype=np.float32) for v in image_volumes), axis=-1)
    fields = np.stack(tuple(np.asarray(v.dataobj, dtype=np.float32) for v in field_volumes), axis=-1)
    memberships = np.asarray(membership_volume.dataobj, dtype=np.float32)
    segmentation = np.asarray(segmentation_volume.dataobj, dtype=np.int16)
    shapes = np.asarray(shape_volume.dataobj, dtype=np.float32)
    assert images.shape == (14, 14, 14, 2)
    assert shapes.shape == (14, 14, 14, 11)
    assert np.all(shapes > 0)
    stats = json.loads((SOURCE_DIR / "statistics_inputs.json").read_text())
    params = json.loads((SOURCE_DIR / "membership_parameters.json").read_text())
    old_centers = np.stack((
        np.asarray(stats["initial_centroids"], dtype=np.float32),
        np.asarray([(11 - k) / 12.0 for k in range(11)], dtype=np.float32),
    ))
    exponents = np.repeat(np.asarray(stats["category_exponents"], dtype=np.float32)[None, :], 2, axis=0)
    field_result = toads_multichannel_exact_field_refit(images, memberships, segmentation, old_centers, fields)

    fields_path = tmp_path / "two_channel_fields.bin"
    shape_path = tmp_path / "transformed_shapes.bin"
    fields_path.write_bytes(np.asarray(field_result.fields, dtype="<f4").tobytes(order="C"))
    shape_path.write_bytes(np.asarray(shapes, dtype="<f4").tobytes(order="C"))
    reference = json.loads(subprocess.run(
        [str(JAVA_HARNESS), str(fields_path), str(shape_path)],
        check=True, capture_output=True, text=True,
    ).stdout)

    harness_source = ROOT / "java_reference/ToadsP18czoNonzeroShapesDynamicMultichannelMembershipReferenceHarness.java"
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["wrapper_source_sha256"] == WRAPPER_SOURCE_SHA256
    assert reference["power_table_source_sha256"] == POWER_TABLE_SHA256
    assert hashlib.sha256(harness_source.read_bytes()).hexdigest() == CASE["harness_sha256"]
    assert reference["reference_scope"].startswith("selected Java formula harness; fixed source-connected synthetic grow/thin state")
    assert reference["source_order"] == [
        "computeExactCentroids", "computeExactVariances", "computeExactImageVariance", "computeExactMemberships"
    ]
    assert reference["shape_sha256"] == _sha(shapes)
    assert reference["shape_sha256"] == CASE["shape_array_sha256"]
    assert reference["shape_nonzero_count"] == int(np.count_nonzero(shapes))
    assert reference["direct_shape_branch_count"] > 0
    assert reference["neighbor_shape_branch_count"] > 0
    assert reference["shape_only_branch_count"] > 0
    assert reference["smoothing_prior_voxel_class_count"] > 0
    assert reference["inverse_shape_prior_voxel_class_count"] > 0

    centers = toads_default_exact_centroid_update(
        images, memberships, old_centers, exponents,
        intensity_field=field_result.fields, smoothness=stats["centroid_smoothness"],
    )
    deviations = toads_default_exact_class_variances(
        images, memberships, centers, segmentation, intensity_field=field_result.fields
    )
    assert _bits(centers[0]) == reference["centroids_channel0_bits"]
    assert _bits(centers[1]) == reference["centroids_channel1_bits"]
    assert _bits(deviations[0]) == reference["deviations_channel0_bits"]
    assert _bits(deviations[1]) == reference["deviations_channel1_bits"]
    image_variances, image_stddev = toads_optimized_image_variance(
        images, memberships[..., 1:], centers[:, 1:], intensity_field=field_result.fields
    )
    assert _bits(image_variances) == reference["image_variance_bits"]
    assert _bits(image_stddev) == reference["image_stddev_bits"]

    categories = np.asarray(params["category_exponents"], dtype=np.float32)
    exponent_by_channel = np.empty((2, 11), dtype=np.float32)
    exponent_by_channel[:, 0] = 1.0
    exponent_by_channel[:, 1:3] = categories[0]
    exponent_by_channel[:, 3:8] = categories[1]
    exponent_by_channel[:, 8:11] = categories[2]
    similarity = toads_optimized_centroid_similarity(centers, deviations, image_variances, params["prior_scale"])
    assert _bits(similarity) == reference["similarity_bits"]
    closest = np.asarray(_load(SOURCE_DIR / "membership_closest_labels.nii.gz").dataobj, dtype=np.int8)
    membership_result = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=memberships,
        centroids=centers,
        exponent_by_channel_class=exponent_by_channel,
        modality_weights=np.asarray(params["modality_weights"], dtype=np.float32),
        class_variance_weights=np.asarray(params["class_variance_weights"], dtype=np.float32),
        image_variances=image_variances,
        transformed_shapes=shapes,
        closest_labels=closest,
        smoothing_priors=np.asarray(params["smoothing_priors"], dtype=np.float32),
        centroid_similarity=similarity,
        atlas_shape_priors=np.asarray(params["atlas_shape_priors"], dtype=np.float32),
        smoothing=params["smoothing"],
        prior_factor=params["prior_factor"],
        prior_coefficient=params["prior_coefficient"],
        intensity_field=field_result.fields,
    )
    assert _sha(membership_result.class_data_energy) == reference["class_energy_sha256"]
    assert _sha(membership_result.memberships) == reference["membership_sha256"]
    assert _sha(membership_result.denominator) == reference["denominator_sha256"]
    assert int(np.asarray(membership_result.objective, dtype=np.float32).view(np.int32)) == reference["objective_bits"]
    assert reference["membership_sha256"] != JAVA_BITS["membership_sha256"]
    assert reference["denominator_sha256"] != JAVA_BITS["denominator_sha256"]
    assert reference["objective_bits"] != JAVA_BITS["objective_bits"]
