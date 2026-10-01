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
from pymacruise.toads_thinning import toads_propagate_growing


ROOT = Path(__file__).parent
CASE_DIR = ROOT / "data/synthetic/toads_p18czo_nonzero_shapes"
CASE = json.loads((CASE_DIR / "case.json").read_text())
BASE_CASE_DIR = ROOT / "data/synthetic/toads_p18czo_multichannel_exact_field"
BASE_CASE = json.loads((BASE_CASE_DIR / "case.json").read_text())
SOURCE_DIR = ROOT / "data/synthetic" / BASE_CASE["source_fixture_directory"]
JAVA_HARNESS = ROOT.parent / "scripts/run_java_toads_p18czo_nonzero_shapes_next_growth_reference.sh"
JAVA_SOURCE_SHA256 = "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
WRAPPER_SOURCE_SHA256 = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
POWER_TABLE_SHA256 = "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).astype(int).tolist()


def _sha(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f4").tobytes(order="C")).hexdigest()


def _decode_labels(encoded: str, shape: tuple[int, int, int], *, skeleton: bool = False) -> np.ndarray:
    values = np.empty(len(encoded), dtype=np.int16)
    for index, value in enumerate(encoded):
        if skeleton and value == ".":
            values[index] = -1
        elif value.isdigit():
            values[index] = int(value)
        else:
            values[index] = ord(value.lower()) - ord("a") + 10
    return values.reshape(shape)


def _encode_labels(values: np.ndarray, *, skeleton: bool = False) -> str:
    output = []
    for value in np.asarray(values).ravel():
        number = int(value)
        output.append("." if skeleton and number < 0 else (str(number) if number < 10 else chr(ord("a") + number - 10)))
    return "".join(output)


def _run_java(fields_path: Path, shapes_path: Path) -> dict:
    return json.loads(subprocess.run(
        [str(JAVA_HARNESS), str(fields_path), str(shapes_path)],
        check=True, capture_output=True, text=True,
    ).stdout)


def test_nonzero_shape_membership_drives_java_referenced_next_growth(tmp_path: Path) -> None:
    channel_names = [entry.split(" / ")[0] for entry in BASE_CASE["channel_order"]]
    field_names = [entry.split(" / ")[1] for entry in BASE_CASE["channel_order"]]
    images_ni = tuple(nib.load(str(SOURCE_DIR / name)) for name in channel_names)
    fields_ni = tuple(nib.load(str(SOURCE_DIR / name)) for name in field_names)
    membership_ni = nib.load(str(SOURCE_DIR / "exact_memberships.nii.gz"))
    segmentation_ni = nib.load(str(BASE_CASE_DIR / BASE_CASE["segmentation"]))
    shapes_ni = nib.load(str(CASE_DIR / CASE["shape_file"]))
    affine = images_ni[0].affine
    for volume in (*images_ni, *fields_ni, membership_ni, segmentation_ni, shapes_ni):
        assert nib.aff2axcodes(volume.affine) == tuple(CASE["axis_codes"])
        np.testing.assert_allclose(volume.affine, affine, atol=0, rtol=0)

    images = np.stack([np.asarray(volume.dataobj, dtype=np.float32) for volume in images_ni], axis=-1)
    fields = np.stack([np.asarray(volume.dataobj, dtype=np.float32) for volume in fields_ni], axis=-1)
    initial_memberships = np.asarray(membership_ni.dataobj, dtype=np.float32)
    segmentation = np.asarray(segmentation_ni.dataobj, dtype=np.int16)
    shapes = np.asarray(shapes_ni.dataobj, dtype=np.float32)
    stats = json.loads((SOURCE_DIR / "statistics_inputs.json").read_text())
    params = json.loads((SOURCE_DIR / "membership_parameters.json").read_text())
    old_centers = np.stack((
        np.asarray(stats["initial_centroids"], dtype=np.float32),
        np.asarray([(11 - k) / 12.0 for k in range(11)], dtype=np.float32),
    ))
    exponents = np.repeat(np.asarray(stats["category_exponents"], dtype=np.float32)[None, :], 2, axis=0)
    field_result = toads_multichannel_exact_field_refit(
        images, initial_memberships, segmentation, old_centers,
        np.stack([np.asarray(volume.dataobj, dtype=np.float32) for volume in fields_ni], axis=-1),
    )
    fields_path = tmp_path / "fitted_fields.bin"
    shapes_path = tmp_path / "nonzero_shapes.bin"
    zero_shapes_path = tmp_path / "zero_shapes.bin"
    fields_path.write_bytes(np.asarray(field_result.fields, dtype="<f4").tobytes(order="C"))
    shapes_path.write_bytes(np.asarray(shapes, dtype="<f4").tobytes(order="C"))
    zero_shapes_path.write_bytes(np.zeros_like(shapes, dtype="<f4").tobytes(order="C"))
    java = _run_java(fields_path, shapes_path)
    zero_java = _run_java(fields_path, zero_shapes_path)

    harness_source = ROOT / "java_reference/ToadsP18czoNonzeroShapesNextGrowthReferenceHarness.java"
    assert hashlib.sha256(harness_source.read_bytes()).hexdigest() == CASE["next_growth_harness_sha256"]
    assert java["source_sha256"] == JAVA_SOURCE_SHA256
    assert java["wrapper_source_sha256"] == WRAPPER_SOURCE_SHA256
    assert java["power_table_source_sha256"] == POWER_TABLE_SHA256
    assert java["source_order"][-4:] == [
        "convergenceCheck.continue", "propagateGrowing.initGrowingLabels", "propagateGrowing.queue", "propagateGrowing.topology"
    ]
    assert java["shape_sha256"] == _sha(shapes) == CASE["shape_array_sha256"]
    assert java["shape_nonzero_count"] == int(np.count_nonzero(shapes))

    centers = toads_default_exact_centroid_update(
        images, initial_memberships, old_centers, exponents,
        intensity_field=field_result.fields, smoothness=stats["centroid_smoothness"],
    )
    deviations = toads_default_exact_class_variances(
        images, initial_memberships, centers, segmentation, intensity_field=field_result.fields,
    )
    image_variances, _ = toads_optimized_image_variance(
        images, initial_memberships[..., 1:], centers[:, 1:], intensity_field=field_result.fields,
    )
    categories = np.asarray(params["category_exponents"], dtype=np.float32)
    exponent_by_channel = np.empty((2, 11), dtype=np.float32)
    exponent_by_channel[:, 0] = 1.0
    exponent_by_channel[:, 1:3] = categories[0]
    exponent_by_channel[:, 3:8] = categories[1]
    exponent_by_channel[:, 8:11] = categories[2]
    similarity = toads_optimized_centroid_similarity(centers, deviations, image_variances, params["prior_scale"])
    closest = np.asarray(nib.load(str(SOURCE_DIR / "membership_closest_labels.nii.gz")).dataobj, dtype=np.int8)
    membership_result = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=initial_memberships,
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
    assert _sha(membership_result.memberships) == java["membership_sha256"]

    pre_shape = tuple(java["shape"][:3])
    seed_segmentation = _decode_labels(java["pre_growth_segmentation"], pre_shape)
    seed_skeleton = _decode_labels(java["pre_growth_skeleton"], pre_shape, skeleton=True)
    seed_ordering = np.asarray(java["pre_growth_ordering_bits"], dtype=np.int32).view(np.float32).reshape(pre_shape)
    seed_volumes = np.asarray(java["pre_growth_volume_bits"], dtype=np.int32).view(np.float32)
    growing = toads_propagate_growing(
        seed_segmentation, seed_skeleton, membership_result.memberships, seed_ordering, seed_volumes,
    )
    events = []
    for event in growing.events:
        inserted = [
            [int(np.asarray(item.value, dtype=np.float32).view(np.int32)), item.x, item.y, item.z, item.label]
            for item in event.inserted_neighbors
        ]
        events.append([
            int(np.asarray(event.popped_item.value, dtype=np.float32).view(np.int32)),
            event.popped_item.x, event.popped_item.y, event.popped_item.z, event.popped_item.label,
            -1 if event.topology_accepted is None else int(event.topology_accepted), inserted,
        ])
    assert events == java["growth_events"]
    assert _encode_labels(growing.segmentation) == java["post_growth_segmentation"]
    assert _encode_labels(growing.skeleton, skeleton=True) == java["post_growth_skeleton"]
    assert _bits(growing.ordering) == java["post_growth_ordering_bits"]
    assert _bits(growing.class_volumes) == java["post_growth_volume_bits"]

    assert java["growth_event_count"] == zero_java["growth_event_count"] == 391
    assert java["membership_sha256"] != zero_java["membership_sha256"]
    assert java["growth_events"] != zero_java["growth_events"]
    assert sum(event[5] == 1 for event in java["growth_events"]) == 50
    assert sum(event[5] == 0 for event in java["growth_events"]) == 299
    assert sum(event[5] == 1 for event in zero_java["growth_events"]) == 50
    assert sum(event[5] == 0 for event in zero_java["growth_events"]) == 323
    assert [event[5] for event in java["growth_events"]] != [event[5] for event in zero_java["growth_events"]]
    # The shape prior changes queue order and which candidates topology rejects;
    # this fixture reaches the same final labels despite the distinct trace.
    assert java["post_growth_segmentation"] == zero_java["post_growth_segmentation"]
