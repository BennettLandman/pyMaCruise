from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.plus_lesion_atlas_init import run_plus_lesion_initial_atlas_alignment_checkpoint
from pymacruise.plus_lesion_crop import run_plus_lesion_crop_checkpoint
import pymacruise.plus_lesion_registered_branch as registered_branch
import pymacruise.plus_lesion_registered_loop as registered_loop
from pymacruise.plus_lesion_registered_recrop import run_plus_lesion_registered_recrop_checkpoint
from pymacruise.plus_lesion_first_iteration import run_plus_lesion_first_iteration_checkpoint
import pymacruise.plus_lesion_first_iteration as first_iteration_stage
from pymacruise.lesion_registered_first_iteration import run_lesion_toads_registered_first_iteration
from pymacruise.lesion_registered_first_iteration_branch import run_lesion_toads_registered_first_iteration_branch
from pymacruise.lesion_atlas_priors import (
    load_lesion_atlas_priors,
    load_lesion_atlas_shape_priors,
    load_lesion_atlas_topology,
)
from pymacruise.lesion_standard13_voxel import STANDARD13_TYPES
from pymacruise.plus_lesion_registered_loop import (
    continue_plus_lesion_registered_loop,
    hydrate_plus_lesion_registered_branch,
    persist_plus_lesion_registered_loop_checkpoint,
)
from pymacruise.plus_lesion_loop_workflow import run_plus_lesion_loop_workflow


ROOT = Path(__file__).parent
SYNTHETIC = ROOT / "data/synthetic/plus_lesion_atlas_init"
ATLAS = SYNTHETIC / "mini-lesion-toads-atlas.txt"


def _run_recrop(tmp_path: Path):
    crop = run_plus_lesion_crop_checkpoint(
        SYNTHETIC / "synthetic_t1_flair.nii.gz",
        SYNTHETIC / "synthetic_topology.nii.gz",
        tmp_path / "crop",
        padding=1,
    )
    alignment = run_plus_lesion_initial_atlas_alignment_checkpoint(
        crop.cropped_path,
        crop.manifest_path,
        ATLAS,
        tmp_path / "alignment",
        init_align_iterations=0,
        pyramid_levels=2,
    )
    return run_plus_lesion_registered_recrop_checkpoint(
        crop.normalized_path,
        crop.manifest_path,
        alignment.manifest_path,
        ATLAS,
        tmp_path / "recrop",
    )


def test_registered_branch_checkpoint_completes_and_resumes(tmp_path: Path, monkeypatch) -> None:
    recrop = _run_recrop(tmp_path)
    captured_initial = {}
    original_initial = first_iteration_stage.run_lesion_toads_registered_initial_membership

    def capture_initial(*args, **kwargs):
        captured_initial["initial"] = original_initial(*args, **kwargs)
        return captured_initial["initial"]

    monkeypatch.setattr(
        first_iteration_stage, "run_lesion_toads_registered_initial_membership", capture_initial
    )
    first_iteration = run_plus_lesion_first_iteration_checkpoint(
        recrop.manifest_path, ATLAS, tmp_path / "first-iteration", iterations=1,
        first_thinning_limit=1000.0, source_zero_variance_nan=False,
    )
    captured = {}
    direct_branch = registered_branch.run_lesion_toads_registered_first_iteration_branch

    def capture_branch(*args, **kwargs):
        captured["branch"] = direct_branch(*args, **kwargs)
        return captured["branch"]

    monkeypatch.setattr(
        registered_branch, "run_lesion_toads_registered_first_iteration_branch", capture_branch
    )
    def unexpected_initial(*args, **kwargs):
        raise AssertionError("P41k must reuse P41i state when a first-iteration manifest is supplied")

    monkeypatch.setattr(
        registered_branch, "run_lesion_toads_registered_initial_membership", unexpected_initial
    )
    result = registered_branch.run_plus_lesion_registered_branch_checkpoint(
        recrop.manifest_path,
        ATLAS,
        tmp_path / "branch",
        iterations=1,
        first_thinning_limit=1000.0,
        source_zero_variance_nan=False,
        first_iteration_manifest_path=first_iteration.manifest_path,
    )
    manifest = json.loads(result.manifest_path.read_text())
    summary = json.loads(result.summary_path.read_text())
    assert manifest["status"] == "complete"
    assert manifest["stage"] == "plus-lesion-registered-branch"
    assert "reuses persisted initial-membership/thinning state" in manifest["equivalence_note"]
    operations = manifest["operation_history"]
    for operation in (
        "statistics_refresh:complete",
        "distance_map_refresh:complete",
        "gm_membership_update:complete",
        "shape_optimizer:complete",
        "shape_registration:complete",
        "growth:complete",
        "relation_refresh:complete",
        "materialize_branch_outputs:complete",
    ):
        assert operation in operations
    assert operations.index("gm_membership_update:complete") < operations.index(
        "shape_registration:running"
    ) < operations.index("convergence_growth_relations:running")
    assert operations.index("shape_optimizer:complete") < operations.index(
        "shape_bounds:0:running"
    ) < operations.index("shape_bounds:12:complete") < operations.index("growth:running")
    assert summary["spatial_convention"] == "RAS"
    assert summary["stop"] is False
    assert summary["growing_events"] > 0
    assert summary["relation_count"] > 0
    stage_order = summary["stage_order"]
    assert stage_order.index("computeVentirclsDistance") < stage_order.index("computeBstemDistance")
    assert stage_order.index("computeBstemDistance") < stage_order.index("computeGMDistance")
    assert summary["stage_order"].index("registerShapes") < summary["stage_order"].index("convergence_check")
    assert len(summary["energy_ratio_state"]) == 7
    assert len(summary["refreshed_centroids"]["class_centroids"]) == 2
    assert len(summary["registration_shape_bounds"]) == 13
    assert {"membership_diagnostics", "distance_state"} <= set(manifest["outputs"])
    modalities = nib.load(str(recrop.modalities_path))
    outputs = [
        nib.load(str(result.memberships_path)),
        nib.load(str(result.segmentation_path)),
        nib.load(str(result.skeleton_path)),
        nib.load(str(result.availability_path)),
        nib.load(str(result.ordering_path)),
        nib.load(str(result.relation_map_path)),
    ]
    assert [image.shape for image in outputs] == [
        (*modalities.shape[:3], 17), modalities.shape[:3], modalities.shape[:3],
        (*modalities.shape[:3], 13), modalities.shape[:3], modalities.shape[:3],
    ]
    assert all(nib.aff2axcodes(image.affine) == ("R", "A", "S") for image in outputs)
    assert all(np.isfinite(np.asarray(image.dataobj)).all() for image in outputs)
    diagnostic = nib.load(manifest["outputs"]["membership_diagnostics"]["path"])
    distance_state = nib.load(manifest["outputs"]["distance_state"]["path"])
    assert diagnostic.shape == (*modalities.shape[:3], 17)
    assert distance_state.shape == (*modalities.shape[:3], 6)
    assert all(nib.aff2axcodes(image.affine) == ("R", "A", "S") for image in (diagnostic, distance_state))
    assert np.isfinite(np.asarray(diagnostic.dataobj)).all()
    result_arrays = [np.asarray(image.dataobj) for image in outputs]
    direct = captured["branch"]
    initial_summary = json.loads(recrop.summary_path.read_text())
    images = np.asarray(nib.load(str(recrop.modalities_path)).dataobj, dtype=np.float32)
    shape_values = np.asarray(nib.load(str(recrop.shapes_path)).dataobj, dtype=np.float32)
    centroid_options = {
        "mode": "none", "modalities": ("T1_MPRAGE", "FLAIR"),
        "lesion_power_exponents": (0.4, 1.2),
        "tissue_power_exponents": np.ones((2, 13), dtype=np.float32),
    }
    transition_direct = run_lesion_toads_registered_first_iteration(
        captured_initial["initial"], images,
        np.asarray(initial_summary["initial_class_centroids"], dtype=np.float32),
        np.asarray(initial_summary["lesion_centroid"], dtype=np.float32),
        np.asarray(initial_summary["black_hole_centroid"], dtype=np.float32),
        shape_values, np.ones(13, dtype=np.float32), np.ones(13, dtype=np.float32),
        centroid_options=centroid_options,
        channel_weight_options={"outlier_distance": 0.1, "norm_type": "Identity"},
        outlier_energy=0.1, lesion_power_exponent=1.2,
        source_zero_variance_nan=False,
        atlas_names=tuple(
            item.name for item in load_lesion_atlas_priors(ATLAS, ("T1", "FLAIR")).structures
        ),
    )
    shape_atlas = load_lesion_atlas_shape_priors(ATLAS)
    atlas_topology = load_lesion_atlas_topology(ATLAS)
    reference_image = nib.load(str(recrop.modalities_path))
    direct_branch_result = run_lesion_toads_registered_first_iteration_branch(
        captured_initial["initial"], transition_direct,
        tuple(shape_atlas.shape_maps[item.name] for item in shape_atlas.structures),
        np.asarray(initial_summary["atlas_transform"], dtype=np.float32),
        image_resolutions_mm=tuple(np.linalg.norm(reference_image.affine[:3, :3], axis=0)),
        template_resolutions_mm=tuple(np.linalg.norm(atlas_topology.image.affine[:3, :3], axis=0)),
        object_types=STANDARD13_TYPES,
        classification=np.asarray(nib.load(str(recrop.topology_path)).dataobj),
        template_labels=np.asarray([item.label for item in shape_atlas.structures], dtype=np.int16),
        lesion_class=10, iteration=1, n_iter_max=2, max_distance=0.0,
        algorithm_mode="diffeomorphism", relation_limit=3, align_shapes=True,
        init_align_iter=1,
    )
    np.testing.assert_array_equal(
        direct.transition.membership.membership.class_memberships,
        transition_direct.membership.membership.class_memberships,
    )
    np.testing.assert_array_equal(
        direct.registration.registration.transform,
        direct_branch_result.registration.registration.transform,
    )
    np.testing.assert_array_equal(
        direct.branch.branch.segmentation,
        direct_branch_result.branch.branch.segmentation,
    )
    refreshed = direct.transition.refreshed_distance_maps
    assert refreshed is not None
    constructor_maps = captured_initial["initial"].distance_maps
    assert any(not np.array_equal(current, prior) for current, prior in (
        (refreshed.ventricle, constructor_maps.ventricle),
        (refreshed.gray_matter, constructor_maps.gray_matter),
        (refreshed.brainstem, constructor_maps.brainstem),
    ))
    np.testing.assert_array_equal(
        np.asarray(distance_state.dataobj),
        np.stack((
            refreshed.ventricle_mask, refreshed.ventricle,
            refreshed.gray_matter_mask, refreshed.gray_matter,
            refreshed.brainstem_mask, refreshed.brainstem,
        ), axis=-1).astype(np.int16),
    )
    np.testing.assert_array_equal(
        result_arrays[0],
        np.concatenate((
            direct.transition.membership.membership.class_memberships,
            direct.transition.membership.membership.lesion_membership[..., None],
            direct.transition.membership.membership.black_hole_membership[..., None],
            direct.transition.membership.membership.outlier_membership[..., None],
            direct.transition.membership.membership.lesion_relation[..., None],
        ), axis=-1),
    )
    np.testing.assert_array_equal(result_arrays[1], direct.branch.branch.segmentation)
    np.testing.assert_array_equal(result_arrays[2], direct.branch.branch.growing.skeleton)
    np.testing.assert_array_equal(result_arrays[3], direct.branch.branch.growing.available.astype(np.uint8))
    np.testing.assert_array_equal(result_arrays[4], direct.branch.branch.growing.ordering)
    np.testing.assert_array_equal(result_arrays[5], direct.branch.branch.relations.relation_map)
    membership = direct.transition.membership.membership
    np.testing.assert_array_equal(
        np.asarray(diagnostic.dataobj),
        np.concatenate((
            membership.denominator[..., None], membership.class_energies,
            membership.lesion_energy[..., None], membership.black_hole_energy[..., None],
            membership.running_objective[..., None],
        ), axis=-1),
    )

    hydrated, registration, distance_state = hydrate_plus_lesion_registered_branch(result.manifest_path)
    np.testing.assert_array_equal(
        hydrated.membership_result.membership.class_memberships,
        membership.class_memberships,
    )
    np.testing.assert_array_equal(hydrated.segmentation, direct.branch.branch.segmentation)
    np.testing.assert_array_equal(hydrated.growing.skeleton, direct.branch.branch.growing.skeleton)
    np.testing.assert_array_equal(hydrated.relations.relation_map, direct.branch.branch.relations.relation_map)
    np.testing.assert_array_equal(registration.transform, direct.registration.registration.transform)
    assert distance_state.shape == (*modalities.shape[:3], 6)

    direct_runner = registered_loop.run_lesion_toads_descriptor_loop
    loop_inputs = {}

    def capture_loop_inputs(branch, initial_registration, **kwargs):
        loop_inputs.update(branch=branch, registration=initial_registration, **kwargs)
        return direct_runner(branch, initial_registration, **kwargs)

    monkeypatch.setattr(registered_loop, "run_lesion_toads_descriptor_loop", capture_loop_inputs)
    loop = continue_plus_lesion_registered_loop(
        result.manifest_path, ATLAS, n_iter_max=2, thinning_limit=1000.0,
        source_zero_variance_nan=False,
    )
    direct_loop = direct_runner(
        direct.branch.branch,
        direct.registration.registration,
        **{key: value for key, value in loop_inputs.items() if key not in {"branch", "registration"}},
    )
    assert len(loop.passes) == 1
    assert loop.passes[0].iteration == 2
    assert loop.passes[0].stage_order
    np.testing.assert_array_equal(loop.final_branch.segmentation, direct_loop.final_branch.segmentation)
    np.testing.assert_array_equal(
        loop.final_branch.membership_result.membership.class_memberships,
        direct_loop.final_branch.membership_result.membership.class_memberships,
    )
    np.testing.assert_array_equal(
        loop.final_branch.registration_handoff.transform,
        direct_loop.final_branch.registration_handoff.transform,
    )

    def unexpected(*args, **kwargs):
        raise AssertionError("valid resume must skip branch computation")

    monkeypatch.setattr(registered_branch, "run_lesion_toads_registered_initial_membership", unexpected)
    resumed = registered_branch.run_plus_lesion_registered_branch_checkpoint(
        recrop.manifest_path,
        ATLAS,
        tmp_path / "branch",
        iterations=1,
        first_thinning_limit=1000.0,
        source_zero_variance_nan=False,
        first_iteration_manifest_path=first_iteration.manifest_path,
        resume=True,
    )
    assert resumed.resumed


def test_registered_branch_rejects_mismatch_and_repairs_corruption(tmp_path: Path) -> None:
    recrop = _run_recrop(tmp_path)
    output = tmp_path / "branch"
    result = registered_branch.run_plus_lesion_registered_branch_checkpoint(
        recrop.manifest_path, ATLAS, output, iterations=1,
        first_thinning_limit=1000.0, source_zero_variance_nan=False,
    )
    with pytest.raises(ValueError, match="settings differ"):
        registered_branch.run_plus_lesion_registered_branch_checkpoint(
            recrop.manifest_path, ATLAS, output, iterations=2,
            first_thinning_limit=1000.0, source_zero_variance_nan=False, resume=True,
        )
    result.segmentation_path.write_bytes(b"corrupt")
    repaired = registered_branch.run_plus_lesion_registered_branch_checkpoint(
        recrop.manifest_path, ATLAS, output, iterations=1,
        first_thinning_limit=1000.0, source_zero_variance_nan=False, resume=True,
    )
    assert not repaired.resumed
    assert nib.load(str(repaired.segmentation_path)).shape == nib.load(str(recrop.modalities_path)).shape[:3]


def test_registered_loop_checkpoint_persists_resumes_and_continues(tmp_path: Path, monkeypatch) -> None:
    recrop = _run_recrop(tmp_path)
    initial = registered_branch.run_plus_lesion_registered_branch_checkpoint(
        recrop.manifest_path, ATLAS, tmp_path / "initial-branch", iterations=1,
        first_thinning_limit=1000.0, source_zero_variance_nan=False,
    )
    captured_loop = {}
    original_continue = registered_loop.continue_plus_lesion_registered_loop

    def capture_continue(*args, **kwargs):
        captured_loop["result"] = original_continue(*args, **kwargs)
        return captured_loop["result"]

    monkeypatch.setattr(registered_loop, "continue_plus_lesion_registered_loop", capture_continue)
    loop = persist_plus_lesion_registered_loop_checkpoint(
        initial.manifest_path, ATLAS, tmp_path / "loop-2", n_iter_max=99,
        thinning_limit=1000.0, max_passes=1, source_zero_variance_nan=False, resume=True,
    )
    manifest = json.loads(loop.manifest_path.read_text())
    summary = json.loads(loop.summary_path.read_text())
    assert manifest["status"] == "complete"
    assert manifest["stage"] == "plus-lesion-registered-loop"
    assert "progress_log" in manifest
    assert summary["completed_iteration"] == 2
    assert summary["stop"] is False
    progress = json.loads(Path(manifest["progress_log"]["path"]).read_text())
    assert progress["run_signature"] == manifest["run_signature"]
    assert len(progress["attempts"]) == 1
    assert progress["attempts"][0]["status"] == "complete"
    assert [event["iteration"] for event in progress["attempts"][0]["events"]] == [2]
    assert progress["attempts"][0]["events"][0]["stage_order"]
    phases = progress["attempts"][0]["phases"]
    pass_phases = [phase for phase in phases if phase["scope"] == "pass"]
    assert [phase["stage"] for phase in pass_phases] == [
        stage
        for stage in (
            "registered_shape_priors",
            "continuing_segmentation_statistics_memberships",
            "registration_and_branching",
        )
        for _ in ("started", "completed")
    ]
    assert [phase["status"] for phase in pass_phases] == ["started", "completed"] * 3
    assert all(
        phase["elapsed_seconds"] >= 0
        for phase in pass_phases
        if "elapsed_seconds" in phase
    )
    transition_phases = [phase for phase in phases if phase["scope"] == "continuing"]
    assert [phase["stage"] for phase in transition_phases] == [
        "propagate_thinning", "propagate_thinning",
        "refresh_statistics", "refresh_statistics",
        "refresh_distance_maps", "refresh_distance_maps",
        "update_memberships", "update_memberships",
    ]
    assert [phase["status"] for phase in transition_phases] == ["started", "completed"] * 4
    assert nib.aff2axcodes(nib.load(str(loop.segmentation_path)).affine) == ("R", "A", "S")
    direct_loop = captured_loop["result"]
    refreshed = direct_loop.final_distance_maps
    assert refreshed is not None
    persisted_distances = nib.load(manifest["outputs"]["distance_state"]["path"])
    assert nib.aff2axcodes(persisted_distances.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(
        np.asarray(persisted_distances.dataobj),
        np.stack((
            refreshed.ventricle_mask, refreshed.ventricle,
            refreshed.gray_matter_mask, refreshed.gray_matter,
            refreshed.brainstem_mask, refreshed.brainstem,
        ), axis=-1).astype(np.int16),
    )
    stage_order = direct_loop.passes[-1].stage_order
    assert stage_order.index("computeVentirclsDistance") < stage_order.index("computeBstemDistance")
    assert stage_order.index("computeBstemDistance") < stage_order.index("computeGMDistance")
    assert stage_order.index("computeGMDistance") < stage_order.index("computeMemberships:GMflag=true")

    def unexpected(*args, **kwargs):
        raise AssertionError("valid loop resume must not recompute")

    monkeypatch.setattr(registered_loop, "continue_plus_lesion_registered_loop", unexpected)
    resumed = persist_plus_lesion_registered_loop_checkpoint(
        initial.manifest_path, ATLAS, tmp_path / "loop-2", n_iter_max=99,
        thinning_limit=1000.0, max_passes=1, source_zero_variance_nan=False, resume=True,
    )
    assert resumed.resumed
    monkeypatch.undo()
    loop.segmentation_path.write_bytes(b"corrupt")
    repaired = persist_plus_lesion_registered_loop_checkpoint(
        initial.manifest_path, ATLAS, tmp_path / "loop-2", n_iter_max=99,
        thinning_limit=1000.0, max_passes=1, source_zero_variance_nan=False, resume=True,
    )
    assert not repaired.resumed
    assert nib.load(str(repaired.segmentation_path)).shape == nib.load(str(recrop.modalities_path)).shape[:3]

    third = persist_plus_lesion_registered_loop_checkpoint(
        loop.manifest_path, ATLAS, tmp_path / "loop-3", n_iter_max=99,
        thinning_limit=1000.0, max_passes=1, source_zero_variance_nan=False,
    )
    third_summary = json.loads(third.summary_path.read_text())
    assert third_summary["completed_iteration"] == 3


def test_registered_loop_checkpoint_carries_adaptive_fields_across_restart(tmp_path: Path) -> None:
    recrop = _run_recrop(tmp_path)
    initial = registered_branch.run_plus_lesion_registered_branch_checkpoint(
        recrop.manifest_path, ATLAS, tmp_path / "initial-branch-fields", iterations=1,
        first_thinning_limit=1000.0, source_zero_variance_nan=False,
    )
    modalities = nib.load(str(recrop.modalities_path))
    field_values = np.ones(modalities.shape, dtype=np.float32)
    field_values[..., 0] = 1.05
    field_path = tmp_path / "initial-correction-fields.nii.gz"
    nib.save(nib.Nifti1Image(field_values, modalities.affine), str(field_path))

    first = persist_plus_lesion_registered_loop_checkpoint(
        initial.manifest_path, ATLAS, tmp_path / "loop-fields-2", n_iter_max=99,
        thinning_limit=1000.0, max_passes=1, source_zero_variance_nan=False,
        correction_fields_path=field_path, correction_polynomial_degree=1,
        initial_correction_degree=0,
    )
    first_manifest = json.loads(first.manifest_path.read_text())
    first_summary = json.loads(first.summary_path.read_text())
    assert first_manifest["inputs"]["correction_fields"]["sha256"]
    assert first_manifest["settings"]["correction_polynomial_degree"] == 1
    assert first_manifest["outputs"]["correction_fields"]["sha256"]
    stored = nib.load(first_manifest["outputs"]["correction_fields"]["path"])
    assert nib.aff2axcodes(stored.affine) == ("R", "A", "S")
    assert stored.shape == modalities.shape
    assert np.isfinite(np.asarray(stored.dataobj)).all()
    assert first_summary["correction_polynomial_degree"] == 1
    assert first_summary["final_correction_degree"] in (0, 1)

    second = persist_plus_lesion_registered_loop_checkpoint(
        first.manifest_path, ATLAS, tmp_path / "loop-fields-3", n_iter_max=99,
        thinning_limit=1000.0, max_passes=1, source_zero_variance_nan=False,
    )
    second_summary = json.loads(second.summary_path.read_text())
    assert second_summary["completed_iteration"] == 3
    assert second_summary["correction_polynomial_degree"] == 1
    second_manifest_fields = json.loads(second.manifest_path.read_text())["outputs"]["correction_fields"]
    assert second_manifest_fields["sha256"]


def test_registered_loop_rejects_uncompressed_correction_field_input(tmp_path: Path) -> None:
    recrop = _run_recrop(tmp_path)
    initial = registered_branch.run_plus_lesion_registered_branch_checkpoint(
        recrop.manifest_path, ATLAS, tmp_path / "initial-branch-uncompressed-fields", iterations=1,
        first_thinning_limit=1000.0, source_zero_variance_nan=False,
    )
    modalities = nib.load(str(recrop.modalities_path))
    uncompressed = tmp_path / "correction-fields.nii"
    nib.save(nib.Nifti1Image(np.ones(modalities.shape, dtype=np.float32), modalities.affine), uncompressed)
    with pytest.raises(ValueError, match=r"correction fields must use compressed NIfTI.*nii.gz"):
        persist_plus_lesion_registered_loop_checkpoint(
            initial.manifest_path, ATLAS, tmp_path / "loop-uncompressed-fields",
            thinning_limit=1000.0, source_zero_variance_nan=False,
            correction_fields_path=uncompressed,
        )


def test_registered_loop_persists_terminal_classification_for_cap_and_convergence(
    tmp_path: Path, monkeypatch
) -> None:
    recrop = _run_recrop(tmp_path)
    initial = registered_branch.run_plus_lesion_registered_branch_checkpoint(
        recrop.manifest_path, ATLAS, tmp_path / "initial-branch", iterations=1,
        first_thinning_limit=1000.0, source_zero_variance_nan=False,
    )
    captured = []
    original_continue = registered_loop.continue_plus_lesion_registered_loop

    def capture_continue(*args, **kwargs):
        result = original_continue(*args, **kwargs)
        captured.append(result)
        return result

    monkeypatch.setattr(registered_loop, "continue_plus_lesion_registered_loop", capture_continue)
    terminal_cases = (
        (
            "iteration-limit",
            {"n_iter_max": 2, "max_passes": 1},
            "iteration_limit",
        ),
        (
            "convergence",
            {"n_iter_max": 99, "max_passes": 1, "convergence_threshold": 1.0e6},
            "convergence",
        ),
    )
    for name, settings, reason in terminal_cases:
        checkpoint = persist_plus_lesion_registered_loop_checkpoint(
            initial.manifest_path, ATLAS, tmp_path / name,
            thinning_limit=1000.0, source_zero_variance_nan=False, **settings,
        )
        direct = captured[-1]
        final = direct.final_branch.final_classification
        assert direct.final_branch.stop is True
        assert final is not None
        manifest = json.loads(checkpoint.manifest_path.read_text())
        summary = json.loads(checkpoint.summary_path.read_text())
        assert summary["stop"] is True
        assert summary["stop_reason"] == reason
        assert manifest["outputs"]["classification"]["sha256"]
        assert manifest["outputs"]["lesion_subtype"]["sha256"]
        assert checkpoint.classification_path is not None
        assert checkpoint.lesion_subtype_path is not None
        classification = nib.load(str(checkpoint.classification_path))
        subtype = nib.load(str(checkpoint.lesion_subtype_path))
        assert classification.get_data_dtype() == np.dtype(np.int8)
        assert subtype.get_data_dtype() == np.dtype(np.int16)
        assert nib.aff2axcodes(classification.affine) == ("R", "A", "S")
        assert nib.aff2axcodes(subtype.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(np.asarray(classification.dataobj), final.classification)
        np.testing.assert_array_equal(np.asarray(subtype.dataobj), final.lesion_classification)
        with pytest.raises(registered_loop.RegisteredBranchStateError, match="terminal LesionTOADS checkpoint"):
            hydrate_plus_lesion_registered_branch(checkpoint.manifest_path)

    def unexpected(*args, **kwargs):
        raise AssertionError("valid terminal resume must reuse its hashed outputs")

    monkeypatch.setattr(registered_loop, "continue_plus_lesion_registered_loop", unexpected)
    resumed = persist_plus_lesion_registered_loop_checkpoint(
        initial.manifest_path, ATLAS, tmp_path / "convergence",
        n_iter_max=99, thinning_limit=1000.0, max_passes=1,
        source_zero_variance_nan=False, convergence_threshold=1.0e6, resume=True,
    )
    assert resumed.resumed
    assert resumed.classification_path is not None


def test_registered_loop_cli_routes_checkpoint_settings(tmp_path: Path, monkeypatch) -> None:
    from pymacruise.cli import main

    captured = {}

    def fake_run(manifest, descriptor, output, **kwargs):
        captured.update(manifest=manifest, descriptor=descriptor, output=output, **kwargs)
        return registered_branch.PlusLesionRegisteredBranchResult(
            tmp_path / "manifest.json",
            *(tmp_path / f"{name}.nii.gz" for name in (
                "memberships", "segmentation", "skeleton", "availability", "ordering", "relation",
            )),
            tmp_path / "summary.json",
            True,
        )

    monkeypatch.setattr(
        registered_loop, "persist_plus_lesion_registered_loop_checkpoint", fake_run
    )
    assert main([
        "plus-lesion-registered-loop",
        "--registered-checkpoint-manifest", "branch.json",
        "--atlas-descriptor", str(ATLAS),
        "--output-dir", "loop-out",
        "--iterations", "20",
        "--thinning-limit", "4",
        "--max-passes", "2",
        "--finite-zero-variance",
        "--correction-fields", "fields.nii.gz",
        "--correction-polynomial-degree", "2",
        "--initial-correction-degree", "1",
        "--resume",
    ]) == 0
    assert captured["n_iter_max"] == 20
    assert captured["thinning_limit"] == 4.0
    assert captured["max_passes"] == 2
    assert captured["source_zero_variance_nan"] is False
    assert captured["resume"] is True
    assert captured["correction_fields_path"] == Path("fields.nii.gz")
    assert captured["correction_polynomial_degree"] == 2
    assert captured["initial_correction_degree"] == 1


def test_source_faithful_zero_variance_nan_continues_through_registered_branch(tmp_path: Path) -> None:
    recrop = _run_recrop(tmp_path)
    output = tmp_path / "source-faithful"
    result = registered_branch.run_plus_lesion_registered_branch_checkpoint(
        recrop.manifest_path, ATLAS, output,
        iterations=1, first_thinning_limit=1000.0,
    )
    manifest = json.loads(result.manifest_path.read_text())
    summary = json.loads(result.summary_path.read_text())
    memberships = nib.load(manifest["outputs"]["memberships"]["path"])
    diagnostics = nib.load(manifest["outputs"]["membership_diagnostics"]["path"])
    assert manifest["status"] == "complete"
    assert summary["source_zero_variance_nan"] is True
    assert summary["transition_objective"] == "NaN"
    assert summary["energy_ratio"] == "NaN"
    assert summary["nonfinite_membership_values"] == 0
    assert summary["nonfinite_membership_energy_values"] > 0
    assert np.isfinite(np.asarray(memberships.dataobj)).all()
    assert np.count_nonzero(~np.isfinite(np.asarray(diagnostics.dataobj))) == summary["nonfinite_membership_energy_values"]

    resumed = registered_branch.run_plus_lesion_registered_branch_checkpoint(
        recrop.manifest_path, ATLAS, output,
        iterations=1, first_thinning_limit=1000.0, resume=True,
    )
    assert resumed.resumed


def test_persisted_loop_workflow_connects_and_resumes_p41i_p41k_p41m(tmp_path: Path) -> None:
    recrop = _run_recrop(tmp_path)
    recrop_manifest = json.loads(recrop.manifest_path.read_text())
    modalities = nib.load(recrop_manifest["outputs"]["modalities"]["path"])
    fields_path = tmp_path / "workflow-correction-fields.nii.gz"
    nib.save(nib.Nifti1Image(np.ones(modalities.shape, dtype=np.float32), modalities.affine), fields_path)
    output = tmp_path / "loop-workflow"
    kwargs = {
        "iterations": 1,
        "first_thinning_limit": 1000.0,
        "source_zero_variance_nan": False,
        "n_iter_max": 99,
        "thinning_limit": 1000.0,
        "max_passes": 1,
        "correction_fields_path": fields_path,
        "correction_polynomial_degree": 1,
        "initial_correction_degree": 0,
    }
    result = run_plus_lesion_loop_workflow(recrop.manifest_path, ATLAS, output, **kwargs)
    manifest = json.loads(result.manifest_path.read_text())
    branch = json.loads(result.registered_branch_manifest.read_text())
    loop = json.loads(result.registered_loop_manifest.read_text())
    assert loop["outputs"]["correction_fields"]["sha256"]
    assert manifest["correction_fields"]["sha256"]
    assert manifest["status"] == "complete"
    assert manifest["workflow_executed"] is True
    assert manifest["source_equivalence_claim"] is False
    assert set(manifest["stages"]) == {"first_iteration", "registered_branch", "registered_loop"}
    assert branch["inputs"]["first_iteration_manifest"]["path"] == str(result.first_iteration_manifest)
    assert "reuses persisted initial-membership/thinning state" in branch["equivalence_note"]
    assert loop["stage"] == "plus-lesion-registered-loop"
    loop_summary = json.loads(Path(loop["outputs"]["summary"]["path"]).read_text())
    assert loop_summary["completed_iteration"] == 2
    for checkpoint in (result.first_iteration_manifest, result.registered_branch_manifest, result.registered_loop_manifest):
        child = json.loads(checkpoint.read_text())
        for name, output_record in child["outputs"].items():
            path = Path(output_record["path"])
            if name == "summary":
                continue
            image = nib.load(str(path))
            assert nib.aff2axcodes(image.affine) == ("R", "A", "S")

    resumed = run_plus_lesion_loop_workflow(
        recrop.manifest_path, ATLAS, output, resume=True, **kwargs,
    )
    assert resumed.resumed
    resumed_manifest = json.loads(resumed.manifest_path.read_text())
    assert all(row["resumed"] for row in resumed_manifest["stages"].values())

    loop_summary_path = Path(loop["outputs"]["summary"]["path"])
    loop_summary_path.write_bytes(b"corrupt")
    repaired = run_plus_lesion_loop_workflow(
        recrop.manifest_path, ATLAS, output, resume=True, **kwargs,
    )
    repaired_manifest = json.loads(repaired.manifest_path.read_text())
    assert repaired_manifest["stages"]["first_iteration"]["resumed"]
    assert repaired_manifest["stages"]["registered_branch"]["resumed"]
    assert not repaired_manifest["stages"]["registered_loop"]["resumed"]
    assert json.loads(loop_summary_path.read_text())["completed_iteration"] == 2

    terminal = run_plus_lesion_loop_workflow(
        recrop.manifest_path,
        ATLAS,
        tmp_path / "terminal-workflow",
        **{**kwargs, "n_iter_max": 2},
    )
    terminal_manifest = json.loads(terminal.manifest_path.read_text())
    assert terminal_manifest["status"] == "complete"
    assert terminal_manifest["terminal"] is True
    assert terminal_manifest["terminal_stop_reason"] == "iteration_limit"
    assert terminal_manifest["result_scope"].endswith("terminal classification and lesion subtype")
    assert terminal.classification_path is not None
    assert terminal.lesion_subtype_path is not None
    assert nib.aff2axcodes(nib.load(str(terminal.classification_path)).affine) == ("R", "A", "S")
    assert nib.aff2axcodes(nib.load(str(terminal.lesion_subtype_path)).affine) == ("R", "A", "S")


def test_plus_lesion_loop_workflow_cli_routes_source_faithful_default(tmp_path: Path, monkeypatch) -> None:
    from types import SimpleNamespace

    from pymacruise.cli import main
    import pymacruise.plus_lesion_loop_workflow as workflow

    captured = {}

    def fake_run(parent, descriptor, output, **kwargs):
        captured.update(parent=parent, descriptor=descriptor, output=output, **kwargs)
        return SimpleNamespace(
            manifest_path=tmp_path / "run.json", resumed=kwargs["resume"],
            classification_path=None, lesion_subtype_path=None,
        )

    monkeypatch.setattr(workflow, "run_plus_lesion_loop_workflow", fake_run)
    assert main([
        "plus-lesion-loop-run",
        "--registered-recrop-manifest", "recrop.json",
        "--atlas-descriptor", "atlas.txt",
        "--output-dir", "out",
        "--iterations", "3",
        "--n-iter-max", "4",
        "--max-passes", "2",
        "--loop-convergence-threshold", "0.9",
        "--correction-fields", "fields.nii.gz",
        "--correction-polynomial-degree", "2",
        "--initial-correction-degree", "1",
        "--resume",
    ]) == 0
    assert captured["source_zero_variance_nan"] is True
    assert captured["iterations"] == 3
    assert captured["n_iter_max"] == 4
    assert captured["max_passes"] == 2
    assert captured["loop_convergence_threshold"] == 0.9
    assert captured["correction_fields_path"] == Path("fields.nii.gz")
    assert captured["correction_polynomial_degree"] == 2
    assert captured["initial_correction_degree"] == 1
    assert captured["resume"] is True


def test_loop_workflow_records_source_faithful_nan_and_completes_one_pass(tmp_path: Path) -> None:
    recrop = _run_recrop(tmp_path)
    output = tmp_path / "source-faithful-workflow"
    result = run_plus_lesion_loop_workflow(
        recrop.manifest_path, ATLAS, output,
        iterations=1, first_thinning_limit=1000.0, max_passes=1,
    )
    manifest = json.loads(result.manifest_path.read_text())
    first_stage = json.loads((output / "first-iteration/lesion-first-iteration.json").read_text())
    branch_summary = json.loads(Path(result.registered_branch_manifest.parent / "registered_branch_summary.json").read_text())
    loop_manifest = json.loads(result.registered_loop_manifest.read_text())
    loop_summary = json.loads(Path(loop_manifest["outputs"]["summary"]["path"]).read_text())
    assert manifest["status"] == "complete"
    assert {name: stage["status"] for name, stage in manifest["stages"].items()} == {
        "first_iteration": "complete", "registered_branch": "complete", "registered_loop": "complete",
    }
    assert first_stage["status"] == "complete"
    assert json.loads(Path(first_stage["outputs"]["summary"]["path"]).read_text())["initial_energy"] == "NaN"
    assert branch_summary["transition_objective"] == "NaN"
    assert branch_summary["energy_ratio"] == "NaN"
    assert loop_summary["transition_objective"] == "NaN"
    assert loop_summary["energy_ratio"] == "NaN"
    assert loop_summary["nonfinite_membership_energy_values"] > 0


def test_registered_branch_cli_routes_finite_mode(tmp_path: Path, monkeypatch) -> None:
    from pymacruise.cli import main

    captured = {}

    def fake_run(manifest, descriptor, output, **kwargs):
        captured.update(manifest=manifest, descriptor=descriptor, output=output, **kwargs)
        return registered_branch.PlusLesionRegisteredBranchResult(
            tmp_path / "manifest.json",
            *(tmp_path / f"{name}.nii.gz" for name in ("memberships", "segmentation", "skeleton", "availability", "ordering", "relation")),
            tmp_path / "summary.json",
            True,
        )

    monkeypatch.setattr(registered_branch, "run_plus_lesion_registered_branch_checkpoint", fake_run)
    assert main([
        "plus-lesion-registered-branch",
        "--registered-recrop-manifest", "recrop.json",
        "--atlas-descriptor", str(ATLAS),
        "--output-dir", "out",
        "--iterations", "1",
        "--first-iteration-manifest", "first.json",
        "--finite-zero-variance",
        "--resume",
    ]) == 0
    assert captured["source_zero_variance_nan"] is False
    assert captured["iterations"] == 1
    assert captured["first_iteration_manifest_path"] == Path("first.json")
    assert captured["resume"] is True
