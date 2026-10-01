from __future__ import annotations

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
from pymacruise.toads_thinning import toads_propagate_thinning
from pymacruise.toads_iteration import toads_stopping_iteration
from test_toads_multichannel_shape_growth_relations import (
    BASE_CASE, BASE_CASE_DIR, CASE, CASE_DIR, SOURCE_DIR, _decode_labels,
    _membership_state, _run_relation,
)
from test_toads_post_thinning_exact_iteration import _field_input_bytes


ROOT = Path(__file__).parent
THINNING_HARNESS = ROOT.parent / "scripts/run_java_toads_post_registration_next_thinning_reference.sh"
FIELD_HARNESS = ROOT.parent / "scripts/run_java_toads_post_thinning_exact_field_reference.sh"
MEMBERSHIP_HARNESS = ROOT.parent / "scripts/run_java_toads_p18czp_actual_next_loop_membership_reference.sh"
TERMINAL_LABELS_HARNESS = ROOT.parent / "scripts/run_java_toads_p18czp_terminal_labels_reference.sh"
CONTROL_REFERENCE = ROOT / "java_reference/ToadsP18cytRegisteredIterationReferenceHarness.java"


def _sha(values: np.ndarray) -> str:
    import hashlib
    return hashlib.sha256(np.asarray(values, dtype="<f4").tobytes(order="C")).hexdigest()


def _labels(values: str, shape: tuple[int, int, int]) -> np.ndarray:
    out = np.empty(len(values), dtype=np.int8)
    for i, value in enumerate(values):
        out[i] = int(value) if value.isdigit() else ord(value.lower()) - ord("a") + 10
    return out.reshape(shape)


def _thinning_input(labels: np.ndarray, memberships: np.ndarray, shapes: np.ndarray, adjacency: np.ndarray) -> bytes:
    return b"".join((
        np.asarray(labels, dtype=np.int8).tobytes(order="C"),
        np.asarray(memberships, dtype=">f4").tobytes(order="C"),
        np.asarray(shapes, dtype=">f4").tobytes(order="C"),
        np.asarray(adjacency, dtype=np.uint8).tobytes(order="C"),
    ))


def _membership_input(labels, images, fields, memberships, old_centers, shapes, closest) -> bytes:
    return b"".join((
        np.asarray(labels, dtype=np.int8).tobytes(order="C"),
        np.asarray(images, dtype="<f4").tobytes(order="C"),
        np.asarray(fields, dtype="<f4").tobytes(order="C"),
        np.asarray(memberships, dtype="<f4").tobytes(order="C"),
        np.asarray(old_centers, dtype="<f4").tobytes(order="C"),
        np.asarray(shapes, dtype="<f4").tobytes(order="C"),
        np.asarray(closest, dtype=np.int8).tobytes(order="C"),
    ))


def _terminal_labels_input(segmentation, template_labels, previous_classification) -> bytes:
    return b"".join((
        np.asarray(segmentation, dtype=np.int8).tobytes(order="C"),
        np.asarray(template_labels, dtype=np.int8).tobytes(order="C"),
        np.asarray(previous_classification, dtype=np.int8).tobytes(order="C"),
    ))


def test_shape_growth_relation_refresh_flows_through_next_thinning_fields_and_membership(tmp_path: Path) -> None:
    """Continue the selected synthetic two-channel, 11-class state across one loop boundary.

    The atlas transform, transformed shape prior and source fixture images stay
    fixed. The test begins from the existing nonzero-shape membership output,
    carries Java-referenced growth and relation/boundary/closest outputs into
    next thinning, then validates the per-channel degree-one field refit and
    explicit-state statistics/membership formula harness. This is a bounded
    selected-source continuation, not full TOADS or registration validation.
    """
    channel_names = [entry.split(" / ")[0] for entry in BASE_CASE["channel_order"]]
    field_names = [entry.split(" / ")[1] for entry in BASE_CASE["channel_order"]]
    image_volumes = tuple(nib.load(str(SOURCE_DIR / name)) for name in channel_names)
    field_volumes = tuple(nib.load(str(SOURCE_DIR / name)) for name in field_names)
    membership_image = nib.load(str(SOURCE_DIR / "exact_memberships.nii.gz"))
    original_segmentation = nib.load(str(BASE_CASE_DIR / BASE_CASE["segmentation"]))
    shape_image = nib.load(str(CASE_DIR / CASE["shape_file"]))
    affine = image_volumes[0].affine
    for volume in (*image_volumes, *field_volumes, membership_image, original_segmentation, shape_image):
        assert nib.aff2axcodes(volume.affine) == tuple(CASE["axis_codes"])
        np.testing.assert_array_equal(volume.affine, affine)
    images = np.stack([np.asarray(v.dataobj, dtype=np.float32) for v in image_volumes], axis=-1)
    fields = np.stack([np.asarray(v.dataobj, dtype=np.float32) for v in field_volumes], axis=-1)
    initial_memberships = np.asarray(membership_image.dataobj, dtype=np.float32)
    initial_seg = np.asarray(original_segmentation.dataobj, dtype=np.int16)
    shapes = np.asarray(shape_image.dataobj, dtype=np.float32)
    stats = json.loads((SOURCE_DIR / "statistics_inputs.json").read_text())
    params = json.loads((SOURCE_DIR / "membership_parameters.json").read_text())
    old_centers = np.stack((
        np.asarray(stats["initial_centroids"], dtype=np.float32),
        np.asarray([(11 - k) / 12.0 for k in range(11)], dtype=np.float32),
    ))
    exponents = np.repeat(np.asarray(stats["category_exponents"], dtype=np.float32)[None, :], 2, axis=0)

    # Reconstruct the current selected loop state before using the post-growth labels.
    first_field = toads_multichannel_exact_field_refit(images, initial_memberships, initial_seg, old_centers, fields)
    current_centers = toads_default_exact_centroid_update(
        images, initial_memberships, old_centers, exponents,
        intensity_field=first_field.fields, smoothness=stats["centroid_smoothness"],
    )
    current_deviations = toads_default_exact_class_variances(
        images, initial_memberships, current_centers, initial_seg, intensity_field=first_field.fields,
    )
    current_image_variances, _ = toads_optimized_image_variance(
        images, initial_memberships[..., 1:], current_centers[:, 1:], intensity_field=first_field.fields,
    )
    categories = np.asarray(params["category_exponents"], dtype=np.float32)
    exponent_by_channel = np.empty((2, 11), dtype=np.float32)
    exponent_by_channel[:, 0] = 1.0
    exponent_by_channel[:, 1:3] = categories[0]
    exponent_by_channel[:, 3:8] = categories[1]
    exponent_by_channel[:, 8:11] = categories[2]
    current_similarity = toads_optimized_centroid_similarity(
        current_centers, current_deviations, current_image_variances, params["prior_scale"]
    )
    initial_closest = np.asarray(nib.load(str(SOURCE_DIR / "membership_closest_labels.nii.gz")).dataobj, dtype=np.int8)
    current_pass = toads_default_atlas_membership_pass(
        image_values=images, initial_memberships=initial_memberships, centroids=current_centers,
        exponent_by_channel_class=exponent_by_channel,
        modality_weights=np.asarray(params["modality_weights"], dtype=np.float32),
        class_variance_weights=np.asarray(params["class_variance_weights"], dtype=np.float32),
        image_variances=current_image_variances, transformed_shapes=shapes, closest_labels=initial_closest,
        smoothing_priors=np.asarray(params["smoothing_priors"], dtype=np.float32),
        centroid_similarity=current_similarity,
        atlas_shape_priors=np.asarray(params["atlas_shape_priors"], dtype=np.float32),
        smoothing=params["smoothing"], prior_factor=params["prior_factor"],
        prior_coefficient=params["prior_coefficient"], intensity_field=first_field.fields,
    )
    memberships, growth_java, grown = _membership_state(tmp_path)
    np.testing.assert_array_equal(memberships, current_pass.memberships)
    np.testing.assert_array_equal(grown, _decode_labels(growth_java["post_growth_segmentation"], grown.shape))

    # The relation refresh consumes the actual continuing growth state and the current memberships.
    prior_java = _run_relation(initial_seg, memberships, np.zeros_like(memberships, dtype=np.int8), tmp_path, "p18czp_prior")
    prior_boundary = np.asarray(prior_java["boundary_distances"], dtype=np.int8).reshape(memberships.shape)
    refresh_java = _run_relation(grown, memberships, prior_boundary, tmp_path, "p18czp_refresh")
    refresh = toads_refresh_relation_state(grown, memberships, prior_boundary, CASE["relation_limit"])
    refreshed_relation_map = np.asarray(refresh_java["relation_map"], dtype=np.int16).reshape(grown.shape)
    refreshed_relation_rows = np.asarray(refresh_java["relation_membership"], dtype=np.int8)
    refreshed_boundary = np.asarray(refresh_java["boundary_distances"], dtype=np.int8).reshape(memberships.shape)
    refreshed_closest = np.asarray(refresh_java["closest_labels"], dtype=np.int8).reshape((*grown.shape, 4))
    np.testing.assert_array_equal(refreshed_relation_map, refresh.relation.relation_map)
    np.testing.assert_array_equal(refreshed_relation_rows.astype(bool), refresh.relation.relation_membership)
    np.testing.assert_array_equal(refreshed_boundary, refresh.boundary_distances)
    np.testing.assert_array_equal(refreshed_closest, refresh.closest_labels)
    adjacency = refresh.relation.class_adjacency

    # Selected wrapper returns to the next loop at propagateThinning; its default shape, volume ratio and limit apply.
    thinning = toads_propagate_thinning(grown, memberships, adjacency, shapes, volume_ratio=.5, limit=3.0)
    thinning_input = tmp_path / "p18czp_thinning_input.bin"
    thinning_input.write_bytes(_thinning_input(grown, memberships, shapes, adjacency))
    java_thin = json.loads(subprocess.run(
        ["sh", str(THINNING_HARNESS), str(thinning_input), *(str(v) for v in grown.shape), "11"],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert java_thin["source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert java_thin["segmentation"] == "".join(
        str(int(v)) if v < 10 else chr(ord("a") + int(v) - 10) for v in thinning.segmentation.ravel()
    )
    assert java_thin["processed"] == [
        [int(np.asarray(i.value, dtype=np.float32).view(np.int32)), i.x, i.y, i.z, i.label]
        for i in thinning.processed_items
    ]

    # Next-loop correction fit uses the labels after thinning, previous current-loop centroids, current maps and fields.
    refit = toads_multichannel_exact_field_refit(
        images, memberships, thinning.segmentation, current_centers, first_field.fields
    )
    for channel in range(2):
        payload = tmp_path / f"p18czp_field_{channel}.bin"
        payload.write_bytes(_field_input_bytes(
            thinning.segmentation, images[..., channel], memberships, first_field.fields[..., channel], current_centers[channel]
        ))
        java_field = json.loads(subprocess.run(
            ["sh", str(FIELD_HARNESS), str(payload), *(str(v) for v in grown.shape), "11"],
            check=True, capture_output=True, text=True,
        ).stdout)
        field_bits = np.asarray(java_field["field_bits"], dtype=np.int32).view(np.float32).reshape(grown.shape)
        np.testing.assert_allclose(refit.fields[..., channel], field_bits, rtol=0, atol=2e-7)

    next_centers = toads_default_exact_centroid_update(
        images, memberships, current_centers, exponents,
        intensity_field=refit.fields, smoothness=stats["centroid_smoothness"],
    )
    next_deviations = toads_default_exact_class_variances(
        images, memberships, next_centers, thinning.segmentation, intensity_field=refit.fields,
    )
    next_image_variances, _ = toads_optimized_image_variance(
        images, memberships[..., 1:], next_centers[:, 1:], intensity_field=refit.fields,
    )
    next_similarity = toads_optimized_centroid_similarity(
        next_centers, next_deviations, next_image_variances, params["prior_scale"]
    )
    refreshed_membership = toads_default_atlas_membership_pass(
        image_values=images, initial_memberships=memberships, centroids=next_centers,
        exponent_by_channel_class=exponent_by_channel,
        modality_weights=np.asarray(params["modality_weights"], dtype=np.float32),
        class_variance_weights=np.asarray(params["class_variance_weights"], dtype=np.float32),
        image_variances=next_image_variances, transformed_shapes=shapes,
        closest_labels=refreshed_closest,
        smoothing_priors=np.asarray(params["smoothing_priors"], dtype=np.float32),
        centroid_similarity=next_similarity,
        atlas_shape_priors=np.asarray(params["atlas_shape_priors"], dtype=np.float32),
        smoothing=params["smoothing"], prior_factor=params["prior_factor"],
        prior_coefficient=params["prior_coefficient"], intensity_field=refit.fields,
    )
    # Selected formula harness consumes this exact continued state, including Java-refreshed closest labels.
    state_path = tmp_path / "p18czp_actual_next_loop.bin"
    state_path.write_bytes(_membership_input(
        thinning.segmentation, images, refit.fields, memberships, current_centers, shapes, refreshed_closest,
    ))
    java_membership = json.loads(subprocess.run(
        ["sh", str(MEMBERSHIP_HARNESS), str(state_path)], check=True, capture_output=True, text=True
    ).stdout)
    assert java_membership["source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert java_membership["wrapper_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert java_membership["power_table_source_sha256"] == "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
    np.testing.assert_array_equal(next_centers.view(np.int32).ravel(), java_membership["centroids_bits"])
    np.testing.assert_array_equal(next_deviations.view(np.int32).ravel(), java_membership["deviations_bits"])
    np.testing.assert_array_equal(next_image_variances.view(np.int32).ravel(), java_membership["image_variance_bits"])
    assert _sha(refreshed_membership.memberships) == java_membership["membership_sha256"]
    assert _sha(refreshed_membership.denominator) == java_membership["denominator_sha256"]
    assert _sha(refreshed_membership.class_data_energy) == java_membership["class_energy_sha256"]
    assert int(np.asarray(refreshed_membership.objective, dtype=np.float32).view(np.int32)) == java_membership["objective_bits"]
    assert np.count_nonzero(refreshed_boundary) > 0
    assert np.count_nonzero(refreshed_closest[..., 1:] >= 0) > 0
    assert java_membership["reference_scope"].startswith("selected Java formula harness on explicitly supplied next-loop state")

    # Terminal control and label transfer are independent of atlas-transform
    # optimization: preserve the actual refreshed membership objective and
    # its prior objective, while treating the registration result as fixed.
    previous_energy = np.float32(current_pass.objective)
    current_energy = np.float32(refreshed_membership.objective)
    ratio = np.float32(2.0) * np.float32(previous_energy - current_energy)
    ratio = np.float32(np.abs(np.float32(ratio / np.float32(previous_energy + current_energy))))
    threshold = np.nextafter(ratio, np.float32(np.inf))
    previous_classification = np.full(thinning.segmentation.shape, -1, dtype=np.int8)
    template_labels = np.asarray([0, 11, 22, 33, 44, 55, 66, 77, 88, 99, 110], dtype=np.int8)
    terminal = toads_stopping_iteration(
        previous_energy, current_energy, 3, 5, threshold,
        thinning.segmentation, thinning.skeleton, refreshed_membership.memberships,
        np.asarray(growth_java["post_growth_ordering_bits"], dtype=np.int32).view(np.float32).reshape(grown.shape),
        thinning.class_volumes, template_labels, previous_classification,
        do_growing=False,
    )
    assert terminal.control.stop and not terminal.control.run_growing
    assert terminal.growing is None

    control_classes = tmp_path / "p18czp_control_classes"
    subprocess.run(["javac", "-d", str(control_classes), str(CONTROL_REFERENCE)], check=True, capture_output=True, text=True)
    control_java = json.loads(subprocess.run(
        ["java", "-cp", str(control_classes), "ToadsP18cytRegisteredIterationReferenceHarness",
         str(int(np.asarray([previous_energy], dtype=np.float32).view(np.int32)[0])),
         str(int(np.asarray([current_energy], dtype=np.float32).view(np.int32)[0])),
         str(int(np.asarray([threshold], dtype=np.float32).view(np.int32)[0])), "3", "5", "false", "true"],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert control_java["source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert control_java["ratio_bits"] == int(np.asarray([terminal.control.energy_ratio]).view(np.int32)[0])
    assert control_java["stop"] is True and control_java["run_growing"] is False
    assert control_java["order"] == [
        "computeExactMemberships", "energyRatio", "atlas.registerShapes",
        "convergenceCheck", "iterationLimitCheck", "updateClassificationFromLabels",
    ]

    terminal_input = tmp_path / "p18czp_terminal_labels.bin"
    terminal_input.write_bytes(_terminal_labels_input(
        thinning.segmentation, template_labels, previous_classification,
    ))
    java_terminal = json.loads(subprocess.run(
        ["sh", str(TERMINAL_LABELS_HARNESS), str(terminal_input), *(str(v) for v in grown.shape), "11"],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert java_terminal["source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    np.testing.assert_array_equal(terminal.classification.ravel(), java_terminal["classification"])
    assert terminal.changed_classification_voxels == java_terminal["changed_voxels"]
