from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.cli import build_parser, main
from pymacruise.image_io import is_axis_aligned_ras_grid
from pymacruise import base_workflow as base_workflow_module
from pymacruise.base_workflow import (
    BaseRunOptions,
    _implementation_signature,
    _materialize_dura_atlas,
    run_base_workflow,
)
from pymacruise import workflow
from pymacruise.surfaces import read_legacy_vtk


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "toads_registered_first_iteration_cycle"
ATLAS = ROOT / "data" / "synthetic" / "toads_registered_field_refit" / "atlas-11.txt"
MASEG = ROOT / "data" / "synthetic" / "base_workflow" / "maseg.nii.gz"
VANDY_ATLAS = ROOT / "data" / "synthetic" / "toads_downsampled_vandy_atlas" / "atlas.txt"


def test_base_cli_regridding_defaults_on_and_can_be_disabled() -> None:
    parser = build_parser()
    required = [
        "base-run", "--t1", "T1.nii.gz", "--brain-mask", "mask.nii.gz",
        "--maseg", "MaSeg.nii.gz", "--atlas-descriptor", "atlas.txt",
        "--output-dir", "run",
    ]
    assert parser.parse_args(required).regrid_inputs is True
    assert parser.parse_args([*required, "--no-regrid-inputs"]).regrid_inputs is False


def test_base_implementation_signature_tracks_python_source_not_notes(tmp_path) -> None:
    source_root = tmp_path / "package"
    source_root.mkdir()
    module = source_root / "module.py"
    notes = source_root / "notes.md"
    module.write_text("VALUE = 1\n", encoding="utf-8")
    notes.write_text("first note\n", encoding="utf-8")
    first = _implementation_signature(source_root)
    assert _implementation_signature(source_root) == first

    notes.write_text("changed note\n", encoding="utf-8")
    assert _implementation_signature(source_root) == first
    module.write_text("VALUE = 2\n", encoding="utf-8")
    assert _implementation_signature(source_root) != first


def test_base_atlas_materializer_accepts_selected_vandy_structure_names(tmp_path) -> None:
    derived = _materialize_dura_atlas(VANDY_ATLAS, tmp_path / "dura-atlas")
    text = derived.read_text()
    assert "Structure: SulcalCSF" in text
    assert "Structure: CerebralGM" in text
    assert "Structure: CerebralWM" in text
    assert (derived.parent / "topology.nii.gz").is_file()
    assert (derived.parent / "shape-10.nii.gz").is_file()


def test_base_config_run_is_end_to_end_restartable_and_stage_checked(
    tmp_path, capsys, monkeypatch
) -> None:
    output = tmp_path / "base-run"
    config_path = tmp_path / "base-workflow.json"
    maseg_values = np.asanyarray(nib.load(str(MASEG)).dataobj).copy()
    # This older workflow fixture only contains brainstem/cerebellum labels.
    # Add a compact MaSeg cortical support so the source-selected base crop and
    # fusion stage can be exercised while retaining all eight cleanup labels.
    maseg_values[4:9, 4:9, 4:9] = 100
    maseg_for_run = tmp_path / "maseg-with-cortex.nii.gz"
    source_t1 = nib.load(str(FIXTURE / "t1.nii.gz"))
    angle = np.deg2rad(8.0)
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle), 0.0],
         [np.sin(angle), np.cos(angle), 0.0],
         [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    oblique_affine = source_t1.affine.copy()
    oblique_affine[:3, :3] = rotation @ oblique_affine[:3, :3]
    t1_for_run = tmp_path / "t1-oblique.nii.gz"
    nib.save(nib.Nifti1Image(np.asarray(source_t1.dataobj), oblique_affine), t1_for_run)
    mask_source = nib.load(str(FIXTURE / "brain-mask.nii.gz"))
    mask_for_run = tmp_path / "brain-mask-oblique.nii.gz"
    nib.save(nib.Nifti1Image(np.asarray(mask_source.dataobj), oblique_affine), mask_for_run)
    nib.save(nib.Nifti1Image(maseg_values.astype(np.int16), oblique_affine), maseg_for_run)
    config_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "profile": "base",
                "t1_path": str(t1_for_run),
                "maseg_path": str(maseg_for_run),
                "brain_mask_path": str(mask_for_run),
                "output_dir": str(output),
                "run_options": {
                    "registration_levels": 1,
                    "registration_iterations": 1,
                    "registration_max_iterations": 1,
                    "toads_iterations": 2,
                    "dura_iterations": 1,
                    "crop_edge_width": 2,
                    "ace_diffusion_iterations": 1,
                    "tgdm_iterations": 1,
                    "tgdm_initial_iterations": 1,
                    "tgdm_inner_iterations": 2,
                    "tgdm_central_iterations": 3,
                    "wm_iso_level": 0.43,
                    "ace_use_updated_csf_in_ace": 0,
                    "source_faithful_refinement": True,
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(workflow, "get_bundled_atlas_descriptor", lambda _name: ATLAS)
    monkeypatch.setattr(workflow, "_ants_available", lambda: True)
    assert main(["workflow-preflight", "--config", str(config_path)]) == 0
    # Simulate an interrupted attempt after preflight and before execution.
    (output / "stages" / "toads-first" / "attempt-01").mkdir(parents=True)
    arguments = ["workflow-run", "--config", str(config_path)]
    assert main(arguments) == 0
    first_summary = capsys.readouterr().out
    assert (output / "preflight.json").is_file()
    assert "Stage check passed: n4" in first_summary
    assert "Stage check passed: toads-first" in first_summary
    assert "Stage check passed: dura" in first_summary
    assert "Stage check passed: toads-second" in first_summary
    assert "Stage check passed: cruise-core" in first_summary
    assert "Stage check passed: maseg-refinement" in first_summary
    assert "Stage check passed: surface-thickness" in first_summary

    manifest_path = output / "run.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["status"] == "complete"
    assert len(manifest["implementation_signature"]) == 64
    assert manifest["spatial_convention"] == "RAS"
    assert manifest["converged"] is False
    assert manifest["result_interpretation"] == "stage-run-complete-but-toads-did-not-meet-energy-threshold"
    assert manifest["source_equivalence_claim"] is False
    regridded = [
        nib.load(item["path"])
        for item in manifest["derived_resources"]
        if "regridded-inputs" in item["path"] and item["path"].endswith(".nii.gz")
    ]
    assert len(regridded) == 3
    assert all(nib.aff2axcodes(image.affine) == ("R", "A", "S") for image in regridded)
    assert all(image.shape == regridded[0].shape for image in regridded)
    assert manifest["stages"]["dura"]["command"][-6:] == [
        "--membership-level", "0.3",
        "--wm-topology-regularization", "0.05",
        "--gm-topology-regularization", "0.02",
    ]
    assert set(manifest["stages"]) == {
        "n4", "toads-first", "dura", "toads-second", "base-tissue-preparation", "cruise-core",
        "maseg-refinement", "surface-thickness",
    }
    assert all(stage["status"] == "complete" for stage in manifest["stages"].values())
    core_command = manifest["stages"]["cruise-core"]["command"]
    for option, expected in (
        ("--tgdm-initial-iterations", "1"),
        ("--tgdm-inner-iterations", "2"),
        ("--tgdm-central-iterations", "3"),
        ("--wm-iso-level", "0.43"),
        ("--use-updated-csf-in-ace", "0"),
    ):
        assert core_command[core_command.index(option) + 1] == expected
    assert "--source-faithful" in manifest["stages"]["maseg-refinement"]["command"]
    assert manifest["stages"]["toads-first"]["attempt"] == 2
    assert all(stage["outputs"] and stage["dependency_signature"] for stage in manifest["stages"].values())
    preparation_stage = manifest["stages"]["base-tissue-preparation"]
    assert len(preparation_stage["outputs"]) == 22  # 20 RAS maps, provenance, stage log
    assert "--remove-cerebellum" not in preparation_stage["command"]
    preparation_provenance_path = Path(next(
        item["path"] for item in preparation_stage["outputs"]
        if Path(item["path"]).name == "provenance.json"
    ))
    preparation_details = json.loads(preparation_provenance_path.read_text())
    assert preparation_details["parameters"]["removed_voxel_counts"] is None


    assert preparation_details["parameters"]["mask_distinction"]["fusion_mask"] == (
        "MaSeg-derived Atlas_BrainTissueMask"
    )
    core_stage = manifest["stages"]["cruise-core"]
    assert "--remove-cerebellum" not in core_stage["command"]
    assert "--base-preparation-provenance" in core_stage["command"]
    assert "--boundary-csf" in core_stage["command"]
    assert any("fused_gm.nii.gz" in value for value in core_stage["command"])
    assert any("boundary_csf.nii.gz" in value for value in core_stage["command"])
    assert any(
        Path(item["path"]).name == "boundary_csf.nii.gz"
        for item in preparation_stage["outputs"]
    )
    toads_log = output / "stages" / "toads-first" / "attempt-02" / "stage.log"
    assert toads_log.is_file()
    toads_log_text = toads_log.read_text()
    assert "TOADS progress +" in toads_log_text
    assert "initial rigid atlas registration complete" in toads_log_text
    assert "registered TOADS segmentation loop complete" in toads_log_text
    assert "first registered cycle complete; iteration=1, energy_ratio=" in toads_log_text
    assert "registered TOADS iteration 2 complete; energy_ratio=" in toads_log_text
    assert "stop=" in toads_log_text
    assert "elapsed_seconds=" in toads_log_text
    for stage_name in ("toads-first", "toads-second"):
        stop = manifest["stages"][stage_name]["toads_stop"]
        assert stop["iterations"] == 2
        assert stop["energy_threshold_reached"] is False
        assert stop["iteration_limit_reached"] is True
        assert stop["terminal_energy_ratio"] > 0

    reference = regridded[0]
    final_dir = Path(manifest["final_outputs"][0]["path"]).parent
    for filename in ("inner.nii.gz", "central.nii.gz", "outer.nii.gz"):
        image = nib.load(str(final_dir / filename))
        values = np.asanyarray(image.dataobj)
        assert image.shape == reference.shape
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert values.dtype == np.float32
        assert np.isfinite(values).all()

    refined = nib.load(str(Path(manifest["stages"]["maseg-refinement"]["outputs"][0]["path"])))
    assert refined.shape == reference.shape
    assert nib.aff2axcodes(refined.affine) == ("R", "A", "S")
    assert np.issubdtype(refined.get_data_dtype(), np.integer)
    assert np.isfinite(np.asanyarray(refined.dataobj)).all()
    core_provenance = Path(next(
        item["path"] for item in manifest["stages"]["cruise-core"]["outputs"]
        if Path(item["path"]).name == "provenance.json"
    ))
    core_details = json.loads(core_provenance.read_text())
    assert core_details["parameters"]["remove_cerebellum"] is False
    assert core_details["parameters"]["removed_voxel_counts"] is None

    surface_stage = manifest["stages"]["surface-thickness"]
    surface_files = {Path(item["path"]).name: Path(item["path"]) for item in surface_stage["outputs"]}
    assert set(surface_files) == {
        "inner_surface.vtk", "central_surface.vtk", "outer_surface.vtk",
        "inner_surface_thickness.vtk", "thickness_levelset.nii.gz",
        "thickness_summary.json", "provenance.json", "stage.log",
    }
    for name in ("inner_surface.vtk", "central_surface.vtk", "outer_surface.vtk"):
        mesh = read_legacy_vtk(surface_files[name])
        assert mesh.coordinate_frame == "RAS"
        assert mesh.units == "mm"
    measured_inner = read_legacy_vtk(surface_files["inner_surface_thickness.vtk"])
    assert "mean_thickness_mm" in measured_inner.point_data
    thickness_image = nib.load(str(surface_files["thickness_levelset.nii.gz"]))
    assert thickness_image.shape == reference.shape
    assert nib.aff2axcodes(thickness_image.affine) == ("R", "A", "S")
    assert np.isfinite(np.asanyarray(thickness_image.dataobj)).all()
    thickness_summary = json.loads(surface_files["thickness_summary.json"].read_text())
    assert thickness_summary["coordinate_frame"] == "RAS"
    assert thickness_summary["levelset_thickness"]["crossing_voxel_count"] > 0
    assert thickness_summary["surface_thickness"]["nonzero_vertex_count"] > 0

    hashes_before = {
        Path(record["path"]).name: record["sha256"]
        for record in manifest["final_outputs"]
        if Path(record["path"]).name.endswith(".nii.gz")
    }
    assert main([*arguments, "--resume"]) == 0
    resumed_summary = capsys.readouterr().out
    assert resumed_summary.count("Resume: verified") == 8
    updated = json.loads(manifest_path.read_text())
    assert {
        Path(record["path"]).name: record["sha256"]
        for record in updated["final_outputs"]
        if Path(record["path"]).name.endswith(".nii.gz")
    } == hashes_before

    # A new package fingerprint invalidates every checkpoint once, then the
    # upgraded manifest resumes normally under that implementation identity.
    monkeypatch.setattr(base_workflow_module, "_implementation_signature", lambda: "f" * 64)
    assert main([*arguments, "--resume"]) == 0
    code_change_summary = capsys.readouterr().out
    assert "Resume: verified" not in code_change_summary
    assert code_change_summary.count("Stage check passed:") == 8
    code_changed = json.loads(manifest_path.read_text())
    assert code_changed["implementation_signature"] == "f" * 64
    assert code_changed["previous_implementation_signature"] == manifest["implementation_signature"]
    assert code_changed["stages"]["n4"]["attempt"] == 2
    assert code_changed["stages"]["toads-first"]["attempt"] == 3
    assert all(
        code_changed["stages"][name]["attempt"] == 2
        for name in ("dura", "toads-second", "cruise-core", "maseg-refinement")
    )

    assert main([*arguments, "--resume"]) == 0
    code_resume_summary = capsys.readouterr().out
    assert code_resume_summary.count("Resume: verified") == 8
    updated = json.loads(manifest_path.read_text())

    # Simulate interruption/corruption at a checkpoint: a damaged CRUISE output
    # forces that stage and its dependent MaSeg refinement into fresh attempts,
    # while all verified upstream stages are reused.
    damaged_inner = Path(next(
        record["path"] for record in updated["stages"]["cruise-core"]["outputs"]
        if Path(record["path"]).name == "inner.nii.gz"
    ))
    damaged_inner.write_bytes(damaged_inner.read_bytes() + b"interrupted-output")
    assert main([*arguments, "--resume"]) == 0
    restart_summary = capsys.readouterr().out
    assert restart_summary.count("Resume: verified") == 5
    assert "Stage check passed: cruise-core" in restart_summary
    assert "Stage check passed: maseg-refinement" in restart_summary
    assert "Stage check passed: surface-thickness" in restart_summary
    restarted = json.loads(manifest_path.read_text())
    assert restarted["stages"]["cruise-core"]["attempt"] == 3
    assert restarted["stages"]["maseg-refinement"]["attempt"] == 3
    assert restarted["stages"]["surface-thickness"]["attempt"] == 3
    assert {
        Path(record["path"]).name: record["sha256"]
        for record in restarted["final_outputs"]
        if Path(record["path"]).name.endswith(".nii.gz")
    } == hashes_before


def test_base_run_uses_orientation_only_ras_grid_when_regridding_is_disabled(tmp_path, monkeypatch) -> None:
    source_t1 = nib.load(str(FIXTURE / "t1.nii.gz"))
    angle = np.deg2rad(8.0)
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle), 0.0],
         [np.sin(angle), np.cos(angle), 0.0],
         [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    affine = source_t1.affine.copy()
    affine[:3, :3] = rotation @ affine[:3, :3]
    t1_path = tmp_path / "oblique-t1.nii.gz"
    mask_path = tmp_path / "oblique-mask.nii.gz"
    maseg_path = tmp_path / "oblique-maseg.nii.gz"
    nib.save(nib.Nifti1Image(np.asarray(source_t1.dataobj), affine), t1_path)
    mask = nib.load(str(FIXTURE / "brain-mask.nii.gz"))
    nib.save(nib.Nifti1Image(np.asarray(mask.dataobj), affine), mask_path)
    maseg = nib.load(str(MASEG))
    nib.save(nib.Nifti1Image(np.asarray(maseg.dataobj), affine), maseg_path)

    def stop_before_n4(command: list[str]) -> None:
        assert command[0] == "n4"
        raise RuntimeError("stopped after no-regrid preflight")

    monkeypatch.setattr(base_workflow_module, "cli_main", stop_before_n4)
    output = tmp_path / "no-regrid"
    with pytest.raises(RuntimeError, match="stopped after no-regrid preflight"):
        run_base_workflow(
            t1_path=t1_path,
            brain_mask_path=mask_path,
            atlas_descriptor=ATLAS,
            maseg_path=maseg_path,
            output_dir=output,
            options=BaseRunOptions(regrid_inputs=False),
        )
    manifest = json.loads((output / "run.json").read_text())
    assert manifest["options"]["regrid_inputs"] is False
    assert not any("regridded-inputs" in item["path"] for item in manifest["derived_resources"])
    oriented = [item for item in manifest["derived_resources"] if "orientation-only-inputs" in item["path"]]
    assert len(oriented) == 3
    oriented_images = {Path(item["path"]).name: nib.load(item["path"]) for item in oriented}
    assert set(oriented_images) == {"T1.nii.gz", "BrainMask.nii.gz", "MaSeg.nii.gz"}
    assert all(is_axis_aligned_ras_grid(image) for image in oriented_images.values())
    expected = nib.as_closest_canonical(nib.load(str(t1_path)))
    np.testing.assert_array_equal(
        np.asarray(oriented_images["T1.nii.gz"].dataobj), np.asarray(expected.dataobj)
    )
    assert manifest["stages"]["n4"]["status"] == "failed"
