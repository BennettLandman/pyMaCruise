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
from pymacruise.plus_lesion_registered_loop import (
    continue_plus_lesion_registered_loop,
    hydrate_plus_lesion_registered_branch,
)


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
    captured = {}
    direct_branch = registered_branch.run_lesion_toads_registered_first_iteration_branch

    def capture_branch(*args, **kwargs):
        captured["branch"] = direct_branch(*args, **kwargs)
        return captured["branch"]

    monkeypatch.setattr(
        registered_branch, "run_lesion_toads_registered_first_iteration_branch", capture_branch
    )
    result = registered_branch.run_plus_lesion_registered_branch_checkpoint(
        recrop.manifest_path,
        ATLAS,
        tmp_path / "branch",
        iterations=1,
        first_thinning_limit=1000.0,
        source_zero_variance_nan=False,
    )
    manifest = json.loads(result.manifest_path.read_text())
    summary = json.loads(result.summary_path.read_text())
    assert manifest["status"] == "complete"
    assert manifest["stage"] == "plus-lesion-registered-branch"
    assert summary["spatial_convention"] == "RAS"
    assert summary["stop"] is False
    assert summary["growing_events"] > 0
    assert summary["relation_count"] > 0
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


def test_source_faithful_zero_variance_failure_is_recorded(tmp_path: Path) -> None:
    recrop = _run_recrop(tmp_path)
    with pytest.raises(ValueError, match="initial membership objective must be finite"):
        registered_branch.run_plus_lesion_registered_branch_checkpoint(
            recrop.manifest_path,
            ATLAS,
            tmp_path / "source-faithful",
            iterations=1,
            first_thinning_limit=1000.0,
        )
    failure = json.loads((tmp_path / "source-faithful/lesion-registered-branch.json").read_text())["failure"]
    assert failure["source_zero_variance_nan"] is True
    assert failure["nonfinite_initial_energy_values"] > 0
    assert failure["zero_class_standard_deviation_channels"] > 0


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
        "--finite-zero-variance",
        "--resume",
    ]) == 0
    assert captured["source_zero_variance_nan"] is False
    assert captured["iterations"] == 1
    assert captured["resume"] is True
