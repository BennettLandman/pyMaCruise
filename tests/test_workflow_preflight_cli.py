from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np
import pytest

from pymacruise import workflow
from pymacruise.cli import main


def _write_input_set(root: Path, *, include_flair: bool = False) -> dict[str, Path]:
    root.mkdir(parents=True, exist_ok=True)
    # The negative x axis creates a physically valid non-RAS input while all
    # three base images remain exactly aligned before normalization.
    affine = np.diag([-1.0, 1.2, 1.5, 1.0])
    affine[0, 3] = 4.0
    shape = (5, 6, 7)
    t1 = np.arange(np.prod(shape), dtype=np.float32).reshape(shape) + 1
    maseg = np.zeros(shape, dtype=np.uint8)
    maseg[1:4, 1:5, 1:6] = 15
    brain_mask = np.zeros(shape, dtype=np.uint8)
    brain_mask[1:4, 1:5, 1:6] = 255
    paths = {
        "t1_path": root / "t1.nii.gz",
        "maseg_path": root / "maseg.nii.gz",
        "brain_mask_path": root / "brain-mask.nii.gz",
    }
    for key, values in (("t1_path", t1), ("maseg_path", maseg), ("brain_mask_path", brain_mask)):
        nib.save(nib.Nifti1Image(values, affine), paths[key])
    if include_flair:
        flair_affine = np.diag([0.8, 0.9, 1.1, 1.0])
        flair = np.ones((6, 7, 8), dtype=np.float32)
        paths["flair_path"] = root / "flair.nii.gz"
        nib.save(nib.Nifti1Image(flair, flair_affine), paths["flair_path"])
    return paths


def _write_config(root: Path, paths: dict[str, Path], *, profile: str = "base") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    config = {
        "schema_version": 1,
        "profile": profile,
        "t1_path": f"../inputs/{paths['t1_path'].name}",
        "maseg_path": f"../inputs/{paths['maseg_path'].name}",
        "brain_mask_path": f"../inputs/{paths['brain_mask_path'].name}",
        "output_dir": "run-output",
    }
    if profile == "plus":
        config["flair_path"] = f"../inputs/{paths['flair_path'].name}"
        config["plus_options"] = {"boundary_overlook_labels": [120, 121]}
    path = root / "workflow.json"
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return path


def test_workflow_preflight_cli_writes_explicit_nonexecuting_manifest(
    tmp_path, monkeypatch, capsys
) -> None:
    inputs = _write_input_set(tmp_path / "inputs")
    config_path = _write_config(tmp_path / "config", inputs)
    # Exercise successful CLI behavior independent of whether the optional
    # ANTsPy extra is installed in the test environment.
    monkeypatch.setattr(workflow, "_ants_available", lambda: True)

    with pytest.warns(UserWarning, match="reorienting the image to RAS") as orientation_warnings:
        assert main(["workflow-preflight", "--config", str(config_path)]) == 0
    captured = capsys.readouterr()
    assert len(orientation_warnings) == 3
    assert "Dry run only" in captured.out
    assert "no N4 correction or segmentation workflow was executed" in captured.out

    manifest_path = tmp_path / "config" / "run-output" / "run.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["workflow_profile"] == "base"
    assert manifest["configuration_schema_version"] == 1
    assert manifest["dry_run"] is True
    assert manifest["stages"]["pipeline"] == "not_started"
    assert manifest["preflight"]["workflow_executed"] is False
    assert manifest["preflight"]["atlas"]["assets_checked"] > 1
    assert manifest["preflight"]["atlas"]["payloads_materialized"] is False
    assert manifest["preflight"]["orientation_conversions"]["t1"] == {
        "from": ["L", "A", "S"],
        "to": ["R", "A", "S"],
        "reoriented": True,
    }
    assert sorted(path.name for path in manifest_path.parent.iterdir()) == ["run.json"]


def test_plus_config_preflight_accepts_separate_flair_grid_and_options(
    tmp_path, monkeypatch
) -> None:
    inputs = _write_input_set(tmp_path / "inputs", include_flair=True)
    config_path = _write_config(tmp_path / "config", inputs, profile="plus")
    monkeypatch.setattr(workflow, "_ants_available", lambda: True)

    config = workflow.load_workflow_config_file(config_path)
    assert config.profile.value == "plus"
    assert config.plus_options is not None
    assert config.plus_options.boundary_overlook_labels == (120, 121)
    manifest_path, manifest = workflow.preflight_workflow(config)
    assert manifest_path.is_file()
    assert set(manifest["inputs"]) == {"t1", "maseg", "brain_mask", "flair"}
    assert manifest["inputs"]["flair"]["output_axis_codes"] == ["R", "A", "S"]
    assert manifest["config"]["plus_options"]["boundary_overlook_labels"] == [120, 121]
    assert manifest["preflight"]["workflow_executed"] is False


def test_workflow_preflight_cli_subprocess_is_dry_run(tmp_path) -> None:
    inputs = _write_input_set(tmp_path / "inputs")
    config_path = _write_config(tmp_path / "config", inputs)
    fake_site = tmp_path / "fake-site"
    fake_site.mkdir()
    # The preflight checks optional-dependency discoverability but never imports
    # or executes ANTs. A marker module keeps this subprocess test core-only.
    (fake_site / "ants.py").write_text("\"\"\"preflight marker only\"\"\"\n")
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(fake_site), environment.get("PYTHONPATH", "")]
    )
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "pymacruise.cli",
            "workflow-preflight",
            "--config",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env=environment,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "Dry run only" in completed.stdout
    assert "reorienting the image to RAS" in completed.stderr
    manifest_path = tmp_path / "config" / "run-output" / "run.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["dry_run"] is True
    assert manifest["stages"]["pipeline"] == "not_started"
    assert not list(manifest_path.parent.glob("*.nii.gz"))


def test_preflight_output_collision_and_overwrite_preserve_other_files(
    tmp_path, monkeypatch
) -> None:
    inputs = _write_input_set(tmp_path / "inputs")
    config_path = _write_config(tmp_path / "config", inputs)
    monkeypatch.setattr(workflow, "_ants_available", lambda: True)
    assert main(["workflow-preflight", "--config", str(config_path)]) == 0
    output_dir = tmp_path / "config" / "run-output"
    sentinel = output_dir / "existing-output.nii.gz"
    sentinel.write_bytes(b"preserved")

    with pytest.raises(FileExistsError, match="use --overwrite"):
        main(["workflow-preflight", "--config", str(config_path)])
    assert main(["workflow-preflight", "--config", str(config_path), "--overwrite"]) == 0
    assert sentinel.read_bytes() == b"preserved"


def test_preflight_reports_actionable_missing_ants_guidance(monkeypatch) -> None:
    config = workflow.WorkflowConfig(
        "base", "t1.nii.gz", "maseg.nii.gz", "mask.nii.gz", "out"
    )
    monkeypatch.setattr(workflow, "_ants_available", lambda: False)
    with pytest.raises(workflow.WorkflowPreflightError, match=r"pymacruise\[ants\]"):
        workflow._check_workflow_dependencies(config)


def test_workflow_config_rejects_unknown_fields_and_non_array_labels(tmp_path) -> None:
    path = tmp_path / "workflow.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "profile": "plus",
                "t1_path": "t1.nii.gz",
                "maseg_path": "maseg.nii.gz",
                "brain_mask_path": "mask.nii.gz",
                "flair_path": "flair.nii.gz",
                "output_dir": "out",
                "unexpected": True,
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Unknown workflow config field"):
        workflow.load_workflow_config_file(path)

    config = json.loads(path.read_text())
    config.pop("unexpected")
    config["plus_options"] = {"boundary_overlook_labels": "120,121"}
    path.write_text(json.dumps(config), encoding="utf-8")
    with pytest.raises(ValueError, match="must be a JSON array"):
        workflow.load_workflow_config_file(path)
