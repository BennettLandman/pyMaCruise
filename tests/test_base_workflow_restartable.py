from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.cli import main
from pymacruise.base_workflow import _materialize_dura_atlas


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "toads_registered_first_iteration_cycle"
ATLAS = ROOT / "data" / "synthetic" / "toads_registered_field_refit" / "atlas-11.txt"
MASEG = ROOT / "data" / "synthetic" / "base_workflow" / "maseg.nii.gz"
VANDY_ATLAS = ROOT / "data" / "synthetic" / "toads_downsampled_vandy_atlas" / "atlas.txt"


def test_base_atlas_materializer_accepts_selected_vandy_structure_names(tmp_path) -> None:
    derived = _materialize_dura_atlas(VANDY_ATLAS, tmp_path / "dura-atlas")
    text = derived.read_text()
    assert "Structure: SulcalCSF" in text
    assert "Structure: CerebralGM" in text
    assert "Structure: CerebralWM" in text
    assert (derived.parent / "topology.nii.gz").is_file()
    assert (derived.parent / "shape-10.nii.gz").is_file()


def test_base_run_is_end_to_end_restartable_and_stage_checked(tmp_path, capsys) -> None:
    output = tmp_path / "base-run"
    # Simulate a previous process that created an attempt directory but died
    # before it could add that attempt to the run manifest.
    (output / "stages" / "toads-first" / "attempt-01").mkdir(parents=True)
    arguments = [
        "base-run",
        "--t1", str(FIXTURE / "t1.nii.gz"),
        "--brain-mask", str(FIXTURE / "brain-mask.nii.gz"),
        "--atlas-descriptor", str(ATLAS),
        "--maseg", str(MASEG),
        "--output-dir", str(output),
        "--registration-levels", "1",
        "--registration-iterations", "1",
        "--registration-max-iterations", "1",
        "--toads-iterations", "2",
        "--dura-iterations", "1",
        "--crop-edge-width", "2",
        "--ace-diffusion-iterations", "1",
        "--tgdm-iterations", "1",
    ]
    assert main(arguments) == 0
    first_summary = capsys.readouterr().out
    assert "Stage check passed: n4" in first_summary
    assert "Stage check passed: toads-first" in first_summary
    assert "Stage check passed: dura" in first_summary
    assert "Stage check passed: toads-second" in first_summary
    assert "Stage check passed: cruise-core" in first_summary
    assert "Stage check passed: maseg-refinement" in first_summary

    manifest_path = output / "run.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["status"] == "complete"
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
