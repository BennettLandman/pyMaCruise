from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

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
        if "atlasem_descriptor_path" in paths:
            config["atlasem_descriptor_path"] = (
                f"../inputs/{paths['atlasem_descriptor_path'].name}"
            )
    path = root / "workflow.json"
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return path


def _write_atlasem_descriptor(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    affine = np.diag([0.8, 0.9, 1.2, 1.0])
    shape = (5, 6, 7)
    nib.save(nib.Nifti1Image(np.zeros(shape, dtype=np.uint8), affine), root / "topology.nii.gz")
    nib.save(nib.Nifti1Image(np.zeros(shape, dtype=np.float32), affine), root / "background.nii.gz")
    nib.save(nib.Nifti1Image(np.ones(shape, dtype=np.float32), affine), root / "cortex.nii.gz")
    descriptor = root / "atlasem.txt"
    descriptor.write_text(
        "\n".join((
            "Structure Atlas File (edit at your own risks)", "", "Structures:\t2",
            "Background\t1\tmask", "Cortex\t2\tclosed", "", "Topology Atlas:",
            "File:\ttopology.nii.gz", "Dimensions:\t5\t6\t7", "Resolutions:\t0.8\t0.9\t1.2",
            "", "Shape Atlas:\t4", "Dimensions:\t5\t6\t7", "Resolutions:\t0.8\t0.9\t1.2",
            "Orientations:\t1\t1\t6\t3", "Structure:\tBackground", "File:\tbackground.nii.gz",
            "Dimensions:\t5\t6\t7", "Resolutions:\t0.8\t0.9\t1.2", "Center:\t0\t0\t0",
            "Structure:\tCortex", "File:\tcortex.nii.gz", "Dimensions:\t5\t6\t7",
            "Resolutions:\t0.8\t0.9\t1.2", "Center:\t0\t0\t0", "",
            "Registered Shapes (0/1: true/false)", "0\t1", "0\t1", "Intensity Atlas:\t12",
            "T1\t0.1\t0.9", "Intensity Variance:", "T1_MPRAGE\t1.5\t2.5",
            "Lesions:", "T1\t0.2\tFLAIR\t0.3", "Black Holes:", "MPRAGE\t0.05",
            "Slope:\t3.5", "Modality weights:", "T1_SPGR\t0.4\t0.6",
            "Optimized factors:", "MPRAGE_RAW\t1.1\t1.2\t1.3\t1.4", "",
        )),
        encoding="utf-8",
    )
    return descriptor


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


def test_plus_config_preflight_tracks_selected_atlasem_assets_and_relative_path(
    tmp_path, monkeypatch
) -> None:
    inputs = _write_input_set(tmp_path / "inputs", include_flair=True)
    inputs["atlasem_descriptor_path"] = _write_atlasem_descriptor(tmp_path / "inputs")
    config_path = _write_config(tmp_path / "config", inputs, profile="plus")
    monkeypatch.setattr(workflow, "_ants_available", lambda: True)

    config = workflow.load_workflow_config_file(config_path)
    assert config.atlasem_descriptor_path == inputs["atlasem_descriptor_path"].resolve()
    manifest_path, manifest = workflow.preflight_workflow(config)
    assert manifest_path.is_file()
    atlasem = manifest["preflight"]["atlasem"]
    assert atlasem["assets_checked"] == 3
    assert atlasem["payloads_materialized"] is False
    assert {Path(item["path"]).name for item in atlasem["assets"]} == {
        "topology.nii.gz", "background.nii.gz", "cortex.nii.gz"
    }
    assert manifest["config"]["atlasem_descriptor_path"] == str(
        inputs["atlasem_descriptor_path"].resolve()
    )
    assert {key for key in manifest["inputs"] if key.startswith("atlasem_")} == {
        "atlasem_descriptor", "atlasem_topology", "atlasem_shape_Background",
        "atlasem_shape_Cortex",
    }


def test_plus_config_preflight_rejects_missing_selected_atlasem_payload(
    tmp_path, monkeypatch
) -> None:
    inputs = _write_input_set(tmp_path / "inputs", include_flair=True)
    descriptor = _write_atlasem_descriptor(tmp_path / "inputs")
    (tmp_path / "inputs" / "cortex.nii.gz").unlink()
    inputs["atlasem_descriptor_path"] = descriptor
    config = workflow.load_workflow_config_file(
        _write_config(tmp_path / "config", inputs, profile="plus")
    )
    monkeypatch.setattr(workflow, "_ants_available", lambda: True)

    with pytest.raises(FileNotFoundError, match="AtlasEM descriptor payloads are missing"):
        workflow.preflight_workflow(config)


def test_plus_config_preflight_rejects_atlasem_payload_on_different_grid(
    tmp_path, monkeypatch
) -> None:
    inputs = _write_input_set(tmp_path / "inputs", include_flair=True)
    descriptor = _write_atlasem_descriptor(tmp_path / "inputs")
    shifted_affine = np.diag([0.8, 0.9, 1.2, 1.0])
    shifted_affine[0, 3] = 0.5
    nib.save(
        nib.Nifti1Image(np.ones((5, 6, 7), dtype=np.float32), shifted_affine),
        tmp_path / "inputs" / "cortex.nii.gz",
    )
    inputs["atlasem_descriptor_path"] = descriptor
    config = workflow.load_workflow_config_file(
        _write_config(tmp_path / "config", inputs, profile="plus")
    )
    monkeypatch.setattr(workflow, "_ants_available", lambda: True)

    with pytest.raises(ValueError, match="grid does not match"):
        workflow.preflight_workflow(config)


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


def test_preflight_overwrite_refuses_symlinked_manifest(tmp_path, monkeypatch) -> None:
    inputs = _write_input_set(tmp_path / "inputs")
    config_path = _write_config(tmp_path / "config", inputs)
    monkeypatch.setattr(workflow, "_ants_available", lambda: True)
    output_dir = tmp_path / "config" / "run-output"
    output_dir.mkdir()
    external_target = tmp_path / "outside.json"
    original = b"preserve this file\n"
    external_target.write_bytes(original)
    (output_dir / "run.json").symlink_to(external_target)

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        with pytest.raises(FileExistsError, match="symlinked preflight manifest"):
            main(
                [
                    "workflow-preflight",
                    "--config",
                    str(config_path),
                    "--overwrite",
                ]
            )

    assert external_target.read_bytes() == original
    assert (output_dir / "run.json").is_symlink()


def test_preflight_reports_actionable_missing_ants_guidance(monkeypatch) -> None:
    config = workflow.WorkflowConfig(
        "base", "t1.nii.gz", "maseg.nii.gz", "mask.nii.gz", "out"
    )
    monkeypatch.setattr(workflow, "_ants_available", lambda: False)
    with pytest.raises(workflow.WorkflowPreflightError, match=r"pymacruise\[ants\]"):
        workflow._check_workflow_dependencies(config)


def test_ants_native_import_error_gets_actionable_guidance(monkeypatch) -> None:
    config = workflow.WorkflowConfig(
        "base", "t1.nii.gz", "maseg.nii.gz", "mask.nii.gz", "out"
    )
    original_find_spec = workflow.importlib.util.find_spec

    def find_spec(name: str):
        if name == "ants":
            return object()
        return original_find_spec(name)

    def import_module(name: str):
        if name == "ants":
            raise OSError("native library unavailable")
        raise AssertionError(f"unexpected import: {name}")

    monkeypatch.setattr(workflow.importlib.util, "find_spec", find_spec)
    monkeypatch.setattr(workflow.importlib, "import_module", import_module)
    with pytest.raises(workflow.WorkflowPreflightError, match=r"pymacruise\[ants\]"):
        workflow._check_workflow_dependencies(config)


def test_missing_bundled_atlas_payload_has_restore_guidance(tmp_path, monkeypatch) -> None:
    descriptor = tmp_path / "toads_atlas.txt"
    descriptor.write_text("File: missing.nii.gz\n", encoding="utf-8")
    monkeypatch.setattr(
        workflow, "get_bundled_atlas_descriptor", lambda _name: descriptor
    )
    with pytest.raises(FileNotFoundError, match="Restore the packaged atlas data"):
        workflow._workflow_atlas_preflight(workflow.WorkflowProfile.BASE)


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


def test_base_config_rejects_atlasem_descriptor(tmp_path) -> None:
    path = tmp_path / "workflow.json"
    path.write_text(json.dumps({
        "schema_version": 1, "profile": "base", "t1_path": "t1.nii.gz",
        "maseg_path": "maseg.nii.gz", "brain_mask_path": "mask.nii.gz",
        "output_dir": "out", "atlasem_descriptor_path": "atlas.txt",
    }), encoding="utf-8")
    with pytest.raises(ValueError, match="only valid|does not accept"):
        workflow.load_workflow_config_file(path)

    plus_config = {
        "schema_version": 1, "profile": "plus", "t1_path": "t1.nii.gz",
        "maseg_path": "maseg.nii.gz", "brain_mask_path": "mask.nii.gz",
        "flair_path": "flair.nii.gz", "output_dir": "out",
        "plus_options": {"boundary_overlook_labels": "120,121"},
    }
    path.write_text(json.dumps(plus_config), encoding="utf-8")
    with pytest.raises(ValueError, match="must be a JSON array"):
        workflow.load_workflow_config_file(path)

    plus_config["plus_options"] = {}
    plus_config["run_options"] = {"plus": {"unexpected": 1}}
    path.write_text(json.dumps(plus_config), encoding="utf-8")
    with pytest.raises(ValueError, match="Unknown run_options.plus field"):
        workflow.load_workflow_config_file(path)


@pytest.mark.parametrize("profile", ["base", "plus"])
def test_workflow_run_cli_dispatches_config_to_restartable_runner(
    tmp_path, monkeypatch, capsys, profile
) -> None:
    inputs = _write_input_set(tmp_path / "inputs", include_flair=profile == "plus")
    config_path = _write_config(tmp_path / "config", inputs, profile=profile)
    descriptor = tmp_path / "atlas.txt"
    descriptor.write_text("synthetic bundled descriptor path\n", encoding="utf-8")
    monkeypatch.setattr(workflow, "get_bundled_atlas_descriptor", lambda _name: descriptor)
    calls = {}
    expected_manifest = tmp_path / "run" / "run.json"
    configured_output = tmp_path / "config" / "run-output"
    configured_output.mkdir(parents=True)
    (configured_output / "run.json").write_text(
        json.dumps(
            {
                "dry_run": True,
                "preflight": {"workflow_executed": False},
            }
        ),
        encoding="utf-8",
    )

    if profile == "base":
        from pymacruise import base_workflow

        def fake_base_runner(**kwargs):
            calls.update(kwargs)
            return expected_manifest

        monkeypatch.setattr(base_workflow, "run_base_workflow", fake_base_runner)
    else:
        from pymacruise import plus_full_workflow

        def fake_plus_runner(*args, **kwargs):
            calls["args"] = args
            calls.update(kwargs)
            return SimpleNamespace(manifest_path=expected_manifest)

        monkeypatch.setattr(plus_full_workflow, "run_plus_full_workflow", fake_plus_runner)

    assert main(["workflow-run", "--config", str(config_path), "--resume"]) == 0
    assert f"Workflow complete: {expected_manifest}" in capsys.readouterr().out
    assert (configured_output / "preflight.json").is_file()
    assert not (configured_output / "run.json").exists()
    assert calls["resume"] is True
    assert calls["atlas_descriptor" if profile == "base" else "atlas_descriptor_path"] == descriptor
    assert calls["maseg_path"] == inputs["maseg_path"]
    if profile == "plus":
        assert calls["pipeline_options"].boundary_overlook_labels == (120, 121)


def test_workflow_run_forwards_configured_atlasem_descriptor(tmp_path, monkeypatch, capsys) -> None:
    inputs = _write_input_set(tmp_path / "inputs", include_flair=True)
    inputs["atlasem_descriptor_path"] = _write_atlasem_descriptor(tmp_path / "inputs")
    config_path = _write_config(tmp_path / "config", inputs, profile="plus")
    descriptor = tmp_path / "lesion-atlas.txt"
    descriptor.write_text("bundled descriptor", encoding="utf-8")
    monkeypatch.setattr(workflow, "get_bundled_atlas_descriptor", lambda _name: descriptor)
    calls = {}
    expected_manifest = tmp_path / "run" / "run.json"
    output_dir = tmp_path / "config" / "run-output"
    output_dir.mkdir(parents=True)
    from pymacruise import plus_full_workflow

    def fake_runner(**kwargs):
        calls.update(kwargs)
        return SimpleNamespace(manifest_path=expected_manifest)

    monkeypatch.setattr(plus_full_workflow, "run_plus_full_workflow", fake_runner)
    assert main(["workflow-run", "--config", str(config_path)]) == 0
    assert "Workflow complete" in capsys.readouterr().out
    assert calls["atlasem_descriptor_path"] == inputs["atlasem_descriptor_path"].resolve()
