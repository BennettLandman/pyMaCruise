from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.cli import main
from pymacruise.base_workflow import _implementation_signature, _materialize_dura_atlas
from pymacruise import workflow


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "toads_registered_first_iteration_cycle"
ATLAS = ROOT / "data" / "synthetic" / "toads_registered_field_refit" / "atlas-11.txt"
MASEG = ROOT / "data" / "synthetic" / "base_workflow" / "maseg.nii.gz"
VANDY_ATLAS = ROOT / "data" / "synthetic" / "toads_downsampled_vandy_atlas" / "atlas.txt"


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
    config_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "profile": "base",
                "t1_path": str(FIXTURE / "t1.nii.gz"),
                "maseg_path": str(MASEG),
                "brain_mask_path": str(FIXTURE / "brain-mask.nii.gz"),
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

    manifest_path = output / "run.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["status"] == "complete"
    assert len(manifest["implementation_signature"]) == 64
    assert manifest["spatial_convention"] == "RAS"
    assert manifest["converged"] is False
    assert manifest["result_interpretation"] == "stage-run-complete-but-toads-did-not-meet-energy-threshold"
    assert manifest["source_equivalence_claim"] is False
    assert set(manifest["stages"]) == {
        "n4", "toads-first", "dura", "toads-second", "cruise-core",
        "maseg-refinement",
    }
    assert all(stage["status"] == "complete" for stage in manifest["stages"].values())
    assert manifest["stages"]["toads-first"]["attempt"] == 2
    assert all(stage["outputs"] and stage["dependency_signature"] for stage in manifest["stages"].values())
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

    reference = nib.load(str(FIXTURE / "t1.nii.gz"))
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
    assert core_details["parameters"]["remove_cerebellum"] is True
    assert core_details["parameters"]["removed_voxel_counts"] == [8, 8, 8]

    hashes_before = {
        Path(record["path"]).name: record["sha256"]
        for record in manifest["final_outputs"]
        if Path(record["path"]).name.endswith(".nii.gz")
    }
    assert main([*arguments, "--resume"]) == 0
    resumed_summary = capsys.readouterr().out
    assert resumed_summary.count("Resume: verified") == 6
    updated = json.loads(manifest_path.read_text())
    assert {
        Path(record["path"]).name: record["sha256"]
        for record in updated["final_outputs"]
        if Path(record["path"]).name.endswith(".nii.gz")
    } == hashes_before

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
    assert restart_summary.count("Resume: verified") == 4
    assert "Stage check passed: cruise-core" in restart_summary
    assert "Stage check passed: maseg-refinement" in restart_summary
    restarted = json.loads(manifest_path.read_text())
    assert restarted["stages"]["cruise-core"]["attempt"] == 2
    assert restarted["stages"]["maseg-refinement"]["attempt"] == 2
    assert {
        Path(record["path"]).name: record["sha256"]
        for record in restarted["final_outputs"]
        if Path(record["path"]).name.endswith(".nii.gz")
    } == hashes_before
