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
from pymacruise.toads_relations import toads_refresh_relation_state


ROOT = Path(__file__).parent
CASE_DIR = ROOT / "data/synthetic/toads_p18czo_nonzero_shapes"
CASE = json.loads((CASE_DIR / "case.json").read_text())
BASE_CASE_DIR = ROOT / "data/synthetic/toads_p18czo_multichannel_exact_field"
BASE_CASE = json.loads((BASE_CASE_DIR / "case.json").read_text())
SOURCE_DIR = ROOT / "data/synthetic" / BASE_CASE["source_fixture_directory"]
GROW_HARNESS = ROOT.parent / "scripts/run_java_toads_p18czo_nonzero_shapes_next_growth_reference.sh"
RELATION_HARNESS = ROOT / "java_reference/run_toads_p18czc_relation_refresh_reference.sh"
SOURCE_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
WRAPPER_SOURCE_SHA256 = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"


def _sha(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f4").tobytes(order="C")).hexdigest()


def _decode_labels(encoded: str, shape: tuple[int, int, int]) -> np.ndarray:
    labels = np.empty(len(encoded), dtype=np.int8)
    for index, value in enumerate(encoded):
        labels[index] = int(value) if value.isdigit() else ord(value.lower()) - ord("a") + 10
    return labels.reshape(shape)


def _run_relation(segmentation: np.ndarray, memberships: np.ndarray, boundary: np.ndarray, tmp_path: Path, tag: str) -> dict:
    seg_path, mem_path, boundary_path = (tmp_path / f"{tag}_{suffix}.bin" for suffix in ("seg", "mem", "boundary"))
    seg_path.write_bytes(np.asarray(segmentation, dtype=np.int8).tobytes(order="C"))
    mem_path.write_bytes(np.asarray(memberships, dtype="<f4").tobytes(order="C"))
    boundary_path.write_bytes(np.asarray(boundary, dtype=np.int8).tobytes(order="C"))
    nx, ny, nz = segmentation.shape
    return json.loads(subprocess.run(
        ["sh", str(RELATION_HARNESS), str(seg_path), str(mem_path), str(boundary_path),
         str(nx), str(ny), str(nz), str(memberships.shape[-1]), str(CASE["relation_limit"])],
        check=True, capture_output=True, text=True,
    ).stdout)


def _membership_state(tmp_path: Path) -> tuple[np.ndarray, dict, object]:
    channel_names = [entry.split(" / ")[0] for entry in BASE_CASE["channel_order"]]
    field_names = [entry.split(" / ")[1] for entry in BASE_CASE["channel_order"]]
    images_ni = tuple(nib.load(str(SOURCE_DIR / name)) for name in channel_names)
    fields_ni = tuple(nib.load(str(SOURCE_DIR / name)) for name in field_names)
    membership_ni = nib.load(str(SOURCE_DIR / "exact_memberships.nii.gz"))
    segmentation_ni = nib.load(str(BASE_CASE_DIR / BASE_CASE["segmentation"]))
    shape_ni = nib.load(str(CASE_DIR / CASE["shape_file"]))
    affine = images_ni[0].affine
    for volume in (*images_ni, *fields_ni, membership_ni, segmentation_ni, shape_ni):
        assert nib.aff2axcodes(volume.affine) == tuple(CASE["axis_codes"])
        np.testing.assert_allclose(volume.affine, affine, atol=0, rtol=0)

    images = np.stack([np.asarray(volume.dataobj, dtype=np.float32) for volume in images_ni], axis=-1)
    fields = np.stack([np.asarray(volume.dataobj, dtype=np.float32) for volume in fields_ni], axis=-1)
    memberships = np.asarray(membership_ni.dataobj, dtype=np.float32)
    segmentation = np.asarray(segmentation_ni.dataobj, dtype=np.int16)
    shapes = np.asarray(shape_ni.dataobj, dtype=np.float32)
    stats = json.loads((SOURCE_DIR / "statistics_inputs.json").read_text())
    params = json.loads((SOURCE_DIR / "membership_parameters.json").read_text())
    old_centers = np.stack((
        np.asarray(stats["initial_centroids"], dtype=np.float32),
        np.asarray([(11 - k) / 12.0 for k in range(11)], dtype=np.float32),
    ))
    exponents = np.repeat(np.asarray(stats["category_exponents"], dtype=np.float32)[None, :], 2, axis=0)
    field_result = toads_multichannel_exact_field_refit(images, memberships, segmentation, old_centers, fields)
    fields_path = tmp_path / "fields.bin"
    fields_path.write_bytes(np.asarray(field_result.fields, dtype="<f4").tobytes(order="C"))
    (tmp_path / "shapes.bin").write_bytes(np.asarray(shapes, dtype="<f4").tobytes(order="C"))
    growth_java = json.loads(subprocess.run(
        [str(GROW_HARNESS), str(fields_path), str(tmp_path / "shapes.bin")],
        check=True, capture_output=True, text=True,
    ).stdout)

    centers = toads_default_exact_centroid_update(
        images, memberships, old_centers, exponents,
        intensity_field=field_result.fields, smoothness=stats["centroid_smoothness"],
    )
    deviations = toads_default_exact_class_variances(
        images, memberships, centers, segmentation, intensity_field=field_result.fields,
    )
    image_variances, _ = toads_optimized_image_variance(
        images, memberships[..., 1:], centers[:, 1:], intensity_field=field_result.fields,
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
    grown = _decode_labels(growth_java["post_growth_segmentation"], tuple(growth_java["shape"][:3]))
    return membership_result.memberships, growth_java, grown


def test_nonzero_shape_multichannel_growth_refreshes_selected_relation_products(tmp_path: Path) -> None:
    memberships, growth_java, grown_segmentation = _membership_state(tmp_path)
    assert _sha(memberships) == growth_java["membership_sha256"]
    assert growth_java["growth_event_count"] == 391
    assert growth_java["post_growth_segmentation"] == growth_java["pre_growth_segmentation"]
    assert hashlib.sha256((ROOT / "java_reference/ToadsP18czoNonzeroShapesNextGrowthReferenceHarness.java").read_bytes()).hexdigest() == CASE["next_growth_harness_sha256"]
    assert hashlib.sha256((ROOT / "java_reference/ToadsP18czcRelationRefreshReferenceHarness.java").read_bytes()).hexdigest() == CASE["relation_refresh_harness_sha256"]

    prior_segmentation = _decode_labels(growth_java["pre_growth_segmentation"], grown_segmentation.shape)
    zero_boundary = np.zeros_like(memberships, dtype=np.int8)
    prior_java = _run_relation(prior_segmentation, memberships, zero_boundary, tmp_path, "prior")
    prior_python = toads_refresh_relation_state(
        prior_segmentation, memberships, zero_boundary, CASE["relation_limit"]
    )
    assert prior_java["source_sha256"] == SOURCE_SHA256
    np.testing.assert_array_equal(
        np.asarray(prior_java["boundary_distances"], dtype=np.int8).reshape(memberships.shape),
        prior_python.boundary_distances,
    )

    refreshed_java = _run_relation(
        grown_segmentation, memberships, prior_python.boundary_distances, tmp_path, "after_growth"
    )
    refreshed = toads_refresh_relation_state(
        grown_segmentation, memberships, prior_python.boundary_distances, CASE["relation_limit"]
    )
    assert refreshed_java["source_sha256"] == SOURCE_SHA256
    assert refreshed_java["relation_diffusion_iterations"] == refreshed.relation.diffusion_iterations
    assert refreshed_java["boundary_diffusion_passes"] == refreshed.boundary_diffusion_passes

    relation_map = np.asarray(refreshed_java["relation_map"], dtype=np.int16).reshape(grown_segmentation.shape)
    relation_table = np.asarray(refreshed_java["relation_membership"], dtype=np.int8)
    boundary_distances = np.asarray(refreshed_java["boundary_distances"], dtype=np.int8).reshape(memberships.shape)
    closest_labels = np.asarray(refreshed_java["closest_labels"], dtype=np.int8).reshape((*grown_segmentation.shape, 4))
    np.testing.assert_array_equal(relation_map, refreshed.relation.relation_map)
    np.testing.assert_array_equal(relation_table.astype(bool), refreshed.relation.relation_membership)
    np.testing.assert_array_equal(boundary_distances, refreshed.boundary_distances)
    np.testing.assert_array_equal(closest_labels, refreshed.closest_labels)
    assert relation_table.shape[0] == refreshed.relation.relation_membership.shape[0] == 178
    assert relation_map.min() >= 0
    assert np.unique(relation_map).size > 11
    assert np.count_nonzero(boundary_distances) > 0
    assert np.count_nonzero(closest_labels[..., 1:] >= 0) > 0
    assert np.unique(relation_map).size == 178
    assert refreshed.relation.diffusion_iterations == 1
    assert refreshed.boundary_diffusion_passes == 1
    assert refreshed.closest_label_passes == 4
