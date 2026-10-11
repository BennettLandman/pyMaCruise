from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES
from pymacruise.toads_initial_thinning_composition import (
    toads_initial_fielded_membership_to_thinning,
)
from pymacruise.toads_thinning import toads_initialize_thinning_queue


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_initial_thinning_startup"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "ToadsInitialThinningStartupReferenceHarness.java"
HARNESS_SHA256 = "6d1b8e91764c7131ccb3683f88e648b81b9b6d9eb491e3145cb16199f80ec631"
REFERENCE_HELPERS = {
    "ToadsHomeomorphicChangeReferenceHarness.java": "b03aada47656df6fc659fc5c66a88430bfceecaf93ae9df086c3dc540496c23c",
    "ToadsRelationListReferenceHarness.java": "38fcd450ae31b751aa9756c16cb835f9cc4ba8351b6d2b392861840538c5cd77",
    "ToadsRelationDiffusionReferenceHarness.java": "7a1ea89d5b3b35e489b40d81cab8a9e385cad850395c531079078944071bc088",
    "ToadsBoundaryDistancesReferenceHarness.java": "86ab1808b8bd2b943ac1821aa0e524d1aee08b42272f2046b57d1c7d796a7353",
    "ToadsClosestLabelsReferenceHarness.java": "7f4c00475807923692c7bdfad436a8042098add39063aa20809ab5767a3e3963",
}


def _sha(values: bytes) -> str:
    return hashlib.sha256(values).hexdigest()


def _float_digest(values: np.ndarray) -> str:
    return _sha(np.asarray(values, dtype=np.float32).view(np.uint32).astype(">u4").tobytes())


def _byte_digest(values: np.ndarray) -> str:
    return _sha(np.asarray(values, dtype=np.int8).tobytes())


def _byte4_digest(values: np.ndarray) -> str:
    return _sha(np.asarray(values, dtype=np.int8).tobytes())


def _int32_digest(values: np.ndarray) -> str:
    return _sha(np.asarray(values, dtype=np.int32).astype(">i4").tobytes())


def _bool_digest(values: np.ndarray) -> str:
    return _sha("".join("1" if value else "0" for value in np.asarray(values, dtype=bool).ravel()).encode("ascii"))


def _available_digest(values: np.ndarray) -> str:
    return _bool_digest(values)


def _queue_digest(items) -> str:
    text = "".join(
        f"{int(np.asarray([item.value], dtype=np.float32).view(np.int32)[0])},"
        f"{item.x},{item.y},{item.z},{item.label}\n"
        for item in items
    )
    return _sha(text.encode("ascii"))


def _load(name: str, *, dtype: np.dtype | type | None = None) -> tuple[np.ndarray, nib.Nifti1Image]:
    image = nib.load(str(ROOT / name))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert image.header.get_xyzt_units()[0] == "mm"
    values = np.asarray(image.dataobj)
    return (values.astype(dtype, copy=False) if dtype is not None else values), image


def test_selected_startup_relation_membership_then_first_thinning_matches_java() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    for filename, digest in REFERENCE_HELPERS.items():
        helper = Path(__file__).parent / "java_reference" / filename
        assert hashlib.sha256(helper.read_bytes()).hexdigest() == digest
    spec = json.loads((ROOT / "spec.json").read_text())
    labels, label_image = _load("segmentation_ras.nii.gz", dtype=np.uint8)
    memberships, membership_image = _load("initial_memberships_ras.nii.gz", dtype=np.float32)
    images, image = _load("images_ras.nii.gz", dtype=np.float32)
    fields, field_image = _load("fields_ras.nii.gz", dtype=np.float32)
    bias, bias_image = _load("psf_bias_ras.nii.gz", dtype=np.float32)
    shapes, shape_image = _load("transformed_shapes_ras.nii.gz", dtype=np.float32)
    for other in (label_image, membership_image, field_image, bias_image, shape_image):
        np.testing.assert_array_equal(other.affine, image.affine)

    assert list(images.shape) == spec["shape_xyz"] + [2]
    assert labels.shape == tuple(spec["shape_xyz"])
    assert memberships.shape == shapes.shape == tuple(spec["shape_xyz"] + [11])
    assert DEFAULT_ATLAS_CLASS_TYPES == tuple(spec["class_order"])
    assert labels.max() == 10 and len(np.unique(labels)) == 11
    assert np.all(shapes > 0.0) and np.ptp(shapes) > 0.3

    params = spec["parameters"]
    fixed_arrays = {
        name: np.asarray(params[name], dtype=np.float32)
        for name in (
            "centroids", "exponent_by_channel_class", "modality_weights",
            "class_variance_weights", "image_variances", "smoothing_priors",
            "centroid_similarity", "atlas_shape_priors",
        )
    }
    fixed_snapshots = {name: values.copy() for name, values in fixed_arrays.items()}
    arrays = {
        "labels": labels,
        "memberships": memberships,
        "images": images,
        "fields": fields,
        "bias": bias,
        "shapes": shapes,
    }
    snapshots = {name: value.copy() for name, value in arrays.items()}
    channel_fields = [fields[..., channel] for channel in range(2)]
    result = toads_initial_fielded_membership_to_thinning(
        image_values=images,
        segmentation=labels,
        initial_memberships=memberships,
        centroids=fixed_arrays["centroids"],
        exponent_by_channel_class=fixed_arrays["exponent_by_channel_class"],
        modality_weights=fixed_arrays["modality_weights"],
        class_variance_weights=fixed_arrays["class_variance_weights"],
        image_variances=fixed_arrays["image_variances"],
        transformed_shapes=shapes,
        smoothing_priors=fixed_arrays["smoothing_priors"],
        centroid_similarity=fixed_arrays["centroid_similarity"],
        atlas_shape_priors=fixed_arrays["atlas_shape_priors"],
        smoothing=spec["membership_smoothing"],
        prior_factor=spec["prior_factor"],
        prior_coefficient=spec["prior_coefficient"],
        channel_fields=channel_fields,
        limit=spec["relation_limit"],
        volume_ratio=spec["thinning_volume_ratio"],
        thinning_limit=spec["thinning_limit"],
        psf_bias=bias,
        # The startup comparison harness explicitly pre-fills skeleton=-1.
        initial_skeleton=np.full(labels.shape, -1, dtype=np.int16),
    )

    assert result.stage_order == tuple(REFERENCE["stage_order"])
    assert result.startup.field_state.use_field
    for channel, field in enumerate(channel_fields):
        assert result.startup.field_state.field_by_channel[channel] is field

    relation = result.startup.relation_state
    assert _int32_digest(relation.relation.relation_map) == REFERENCE["relation_map_sha256"]
    assert _bool_digest(relation.relation.class_adjacency) == REFERENCE["relation_adjacency_sha256"]
    assert _byte4_digest(relation.boundary_distances) == REFERENCE["boundary_distances_sha256"]
    assert _byte4_digest(relation.closest_labels) == REFERENCE["closest_labels_sha256"]
    assert _float_digest(memberships) == REFERENCE["starting_membership_sha256"]
    assert _float_digest(shapes) == REFERENCE["transformed_shape_samples_sha256"]
    assert REFERENCE["fixed_upstream"].startswith("class centroids, class variance weights")
    assert REFERENCE["registration_transform"].startswith("not computed or consumed")

    membership = result.membership_pass
    assert _float_digest(membership.class_data_energy) == REFERENCE["class_energy_sha256"]
    assert _float_digest(membership.memberships) == REFERENCE["membership_sha256"]
    assert _float_digest(membership.denominator) == REFERENCE["denominator_sha256"]
    assert int(np.asarray(membership.objective, dtype=np.float32).view(np.int32)) == REFERENCE["objective_bits"]

    initial_thinning = toads_initialize_thinning_queue(labels, membership.memberships)
    thinning = result.thinning
    assert initial_thinning.boundary_counts.tolist() == REFERENCE["boundary_counts"]
    np.testing.assert_array_equal(
        initial_thinning.class_volumes.view(np.int32),
        np.asarray(REFERENCE["initial_volume_bits"], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        thinning.minimum_class_volumes.view(np.int32),
        np.asarray(REFERENCE["minimum_volume_bits"], dtype=np.int32),
    )
    assert len(thinning.processed_items) == REFERENCE["processed_count"]
    assert _queue_digest(thinning.processed_items) == REFERENCE["processed_sha256"]
    np.testing.assert_array_equal(
        thinning.class_volumes.view(np.int32),
        np.asarray(REFERENCE["final_volume_bits"], dtype=np.int32),
    )
    assert _byte_digest(thinning.segmentation) == REFERENCE["segmentation_sha256"]
    assert _byte_digest(thinning.skeleton) == REFERENCE["skeleton_sha256"]
    # The selected Java startup harness pre-fills every voxel with -1; only
    # the interior is assigned during initThinningLabels.
    assert np.all(thinning.skeleton[0, :, :] == -1)
    assert np.all(thinning.skeleton[-1, :, :] == -1)
    assert np.all(thinning.skeleton[:, 0, :] == -1)
    assert np.all(thinning.skeleton[:, -1, :] == -1)
    assert np.all(thinning.skeleton[:, :, 0] == -1)
    assert np.all(thinning.skeleton[:, :, -1] == -1)
    assert _available_digest(thinning.available) == REFERENCE["available_sha256"]
    assert len(thinning.queue) == REFERENCE["queue_after_count"]
    queue_after = []
    while thinning.queue:
        item = thinning.queue.pop_first()
        queue_after.append(item)
    assert _queue_digest(queue_after) == REFERENCE["queue_after_sha256"]

    # At least one class crosses the absolute 100-voxel floor and loses volume;
    # this confirms the synthetic case executes eligible source thinning work.
    assert np.any(initial_thinning.class_volumes > thinning.minimum_class_volumes)
    assert np.any(thinning.class_volumes < initial_thinning.class_volumes)
    for name, before in snapshots.items():
        np.testing.assert_array_equal(arrays[name], before)
    for name, before in fixed_snapshots.items():
        np.testing.assert_array_equal(fixed_arrays[name], before)


def test_startup_thinning_rejects_zero_or_misaligned_shape_samples() -> None:
    spec = json.loads((ROOT / "spec.json").read_text())
    labels, _ = _load("segmentation_ras.nii.gz", dtype=np.uint8)
    memberships, _ = _load("initial_memberships_ras.nii.gz", dtype=np.float32)
    images, _ = _load("images_ras.nii.gz", dtype=np.float32)
    fields, _ = _load("fields_ras.nii.gz", dtype=np.float32)
    shapes, _ = _load("transformed_shapes_ras.nii.gz", dtype=np.float32)
    p = spec["parameters"]
    args = dict(
        image_values=images, segmentation=labels, initial_memberships=memberships,
        centroids=np.asarray(p["centroids"], dtype=np.float32),
        exponent_by_channel_class=np.asarray(p["exponent_by_channel_class"], dtype=np.float32),
        modality_weights=np.asarray(p["modality_weights"], dtype=np.float32),
        class_variance_weights=np.asarray(p["class_variance_weights"], dtype=np.float32),
        image_variances=np.asarray(p["image_variances"], dtype=np.float32),
        smoothing_priors=np.asarray(p["smoothing_priors"], dtype=np.float32),
        centroid_similarity=np.asarray(p["centroid_similarity"], dtype=np.float32),
        atlas_shape_priors=np.asarray(p["atlas_shape_priors"], dtype=np.float32),
        smoothing=spec["membership_smoothing"], prior_factor=spec["prior_factor"],
        prior_coefficient=spec["prior_coefficient"],
        channel_fields=[fields[..., 0], fields[..., 1]], limit=spec["relation_limit"],
    )
    with np.testing.assert_raises_regex(ValueError, "transformed_shapes must be finite"):
        toads_initial_fielded_membership_to_thinning(
            **args, transformed_shapes=np.full_like(shapes, np.nan)
        )
    with np.testing.assert_raises_regex(ValueError, "transformed_shapes must be finite and match"):
        toads_initial_fielded_membership_to_thinning(
            **args, transformed_shapes=shapes[..., :10]
        )
