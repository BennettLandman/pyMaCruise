from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise import base_resume_workflow as resume_module
from pymacruise.base_resume_workflow import run_base_resume_from_toads
from pymacruise.cli import main
from pymacruise.surfaces import SurfaceMesh, write_legacy_vtk


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(path: Path, values: np.ndarray, affine: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(values, affine), path)


def _source_run(root: Path, *, shape=(8, 9, 10)) -> tuple[Path, Path]:
    affine = np.diag([1.2, 1.1, 1.4, 1.0])
    t1 = root / "source-t1.nii.gz"
    maseg = root / "maseg.nii.gz"
    brain_mask = root / "brain-mask.nii.gz"
    atlas_descriptor = root / "atlas.txt"
    atlas_asset = root / "atlas-shape.nii.gz"
    base = np.arange(np.prod(shape), dtype=np.float32).reshape(shape)
    _save(t1, base, affine)
    _save(maseg, (base.astype(np.int16) % 5), affine)
    _save(brain_mask, np.ones(shape, dtype=np.uint8), affine)
    atlas_descriptor.write_text("File: atlas-shape.nii.gz\n")
    _save(atlas_asset, np.ones(shape, dtype=np.float32), affine)
    materialized_atlas = root / "resources" / "dura-atlas" / "atlas.txt"
    materialized_atlas.parent.mkdir(parents=True)
    materialized_atlas.write_text("File: atlas-shape.nii.gz\n")
    (materialized_atlas.parent / "atlas-shape.nii.gz").write_bytes(atlas_asset.read_bytes())
    stage_dir = root / "source-stage" / "attempt-01"
    stage_dir.mkdir(parents=True)
    paths = {
        "classification.nii.gz": np.ones(shape, dtype=np.uint16),
        "stripped_t1_gm.nii.gz": np.full(shape, 0.7, dtype=np.float32),
        "stripped_t1_wmfill.nii.gz": np.full(shape, 0.8, dtype=np.float32),
        "stripped_t1_wmmask.nii.gz": np.ones(shape, dtype=np.uint8),
        "stripped_t1_csf.nii.gz": np.full(shape, 0.1, dtype=np.float32),
    }
    records = []
    for name, values in paths.items():
        path = stage_dir / name
        _save(path, values, affine)
        image = nib.load(str(path))
        records.append({
            "path": str(path), "sha256": _hash(path), "bytes": path.stat().st_size,
            "shape": list(image.shape), "axis_codes": ["R", "A", "S"],
        })
    dura_dir = root / "source-dura" / "attempt-01"
    dura_records = []
    for name, source in (("stripped_t1.nii.gz", t1), ("brain_mask.nii.gz", brain_mask)):
        path = dura_dir / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(source.read_bytes())
        image = nib.load(str(path))
        dura_records.append({
            "path": str(path), "sha256": _hash(path), "bytes": path.stat().st_size,
            "shape": list(image.shape), "axis_codes": ["R", "A", "S"],
        })
    prior_outputs = {}
    for stage_name, filename, array in (
        ("n4", "t1-n4.nii.gz", base),
        ("toads-first", "classification.nii.gz", np.ones(shape, dtype=np.uint16)),
    ):
        path = root / f"source-{stage_name}" / "attempt-01" / filename
        _save(path, np.asarray(array), affine)
        prior_outputs[stage_name] = [{"path": str(path), "sha256": _hash(path), "bytes": path.stat().st_size}]
    log = stage_dir / "stage.log"
    log.write_text("completed\n")
    records.append({"path": str(log), "sha256": _hash(log), "bytes": log.stat().st_size})
    manifest = {
        "schema_version": 1, "status": "complete", "spatial_convention": "RAS",
        "run_signature": "a" * 64, "implementation_signature": "b" * 64,
        "workflow_scope": "synthetic/native-grid base route; selected V3 run", "options": {
            "registration_levels": 2, "registration_iterations": 10,
            "post_refresh_registration_iterations": 3, "registration_max_iterations": 4,
            "toads_iterations": 50, "crop_edge_width": 2,
            "dura_iterations": 1, "dura_membership_level": 0.3,
            "ace_diffusion_iterations": 1, "gvf_iterations": 200, "tgdm_iterations": 1,
        },
        "inputs": {
            "t1": {"path": str(t1), "sha256": _hash(t1)},
            "brain_mask": {"path": str(brain_mask), "sha256": _hash(brain_mask)},
            "atlas_descriptor": {"path": str(atlas_descriptor), "sha256": _hash(atlas_descriptor)},
            "maseg": {"path": str(maseg), "sha256": _hash(maseg)},
        },
        "atlas_assets": [{"path": str(atlas_asset), "sha256": _hash(atlas_asset)}],
        "derived_resources": [{"path": str(materialized_atlas), "sha256": _hash(materialized_atlas)}],
        "stages": {
            "n4": {"status": "complete", "dependency_signature": "e" * 64, "outputs": prior_outputs["n4"]},
            "toads-first": {"status": "complete", "dependency_signature": "f" * 64, "outputs": prior_outputs["toads-first"]},
            "dura": {"status": "complete", "dependency_signature": "1" * 64, "outputs": dura_records,
                     "command": ["dura", "--output-dir", str(dura_dir), "--max-iterations", "1",
                                 "--membership-level", "0.3"]},
            "toads-second": {
            "status": "complete", "attempt": 1, "dependency_signature": "c" * 64,
            "command": ["toads-loop", "--t1-n4", str(dura_dir / "stripped_t1.nii.gz"), "--brain-mask", str(dura_dir / "brain_mask.nii.gz"),
                        "--output", "classification.nii.gz", "--dura-output-dir", "dura-inputs",
                        "--atlas-descriptor", str(materialized_atlas), "--registration-levels", "2",
                        "--registration-iterations", "10", "--post-refresh-registration-iterations", "3",
                        "--registration-max-iterations", "4",
                        "--maximum-iterations", "50"],
            "toads_stop": {"iterations": 23, "terminal_energy_ratio": 9e-5,
                            "energy_threshold_reached": True, "iteration_limit_reached": False},
            "converged": True, "outputs": records,
        }},
    }
    source_run = root / "source-workflow" / "run.json"
    source_run.parent.mkdir(parents=True, exist_ok=True)
    source_run.write_text(json.dumps(manifest))
    return source_run, maseg


def _install_fake_cli(monkeypatch, t1_path: Path) -> list[list[str]]:
    calls: list[list[str]] = []
    reference = nib.load(str(t1_path))
    affine = reference.affine
    shape = reference.shape
    crop_shape = tuple(size - 2 for size in shape)

    def fake_cli(command: list[str]) -> int:
        calls.append(command)
        stage = command[0]
        if stage == "segmentation-refine":
            output_dir = Path(command[command.index("--output") + 1]).parent
        else:
            output_dir = Path(command[command.index("--output-dir") + 1])
        output_dir.mkdir(parents=True, exist_ok=True)
        if stage == "base-tissue-preparation":
            full_names = (
                "atlas_gm.nii.gz", "atlas_wm.nii.gz", "atlas_mask.nii.gz", "atlas_brain_tissue_mask.nii.gz",
                "toads_gm_transferred.nii.gz", "toads_wmfill_transferred.nii.gz", "toads_mask_transferred.nii.gz",
            )
            crop_names = (
                "atlas_gm_cropped.nii.gz", "atlas_wm_cropped.nii.gz", "atlas_mask_cropped.nii.gz",
                "atlas_brain_tissue_mask_cropped.nii.gz", "toads_gm_cropped.nii.gz", "toads_wmfill_cropped.nii.gz",
                "toads_mask_cropped.nii.gz", "fused_gm.nii.gz", "fused_wmfill.nii.gz", "fused_mask.nii.gz",
            )
            for name in full_names:
                _save(output_dir / name, np.ones(shape, dtype=np.float32), affine)
            crop_affine = affine.copy()
            crop_affine[:3, 3] += affine[:3, :3] @ np.ones(3)
            for name in crop_names:
                _save(output_dir / name, np.ones(crop_shape, dtype=np.float32), crop_affine)
            (output_dir / "provenance.json").write_text('{"fixture":true}\n')
        elif stage == "toads-cruise-core":
            crop_affine = affine.copy()
            crop_affine[:3, 3] += affine[:3, :3] @ np.ones(3)
            crop_names = ("gm_cropped.nii.gz", "wmfill_cropped.nii.gz", "wmmask_cropped.nii.gz",
                          "wmfill_topology_corrected.nii.gz", "ace_enhanced_gm.nii.gz",
                          "ace_skeleton.nii.gz", "ace_thinned_skeleton.nii.gz")
            for name in crop_names:
                _save(output_dir / name, np.ones(crop_shape, dtype=np.float32), crop_affine)
            _save(output_dir / "gvf.nii.gz", np.ones((*crop_shape, 3), dtype=np.float32), crop_affine)
            for name in ("inner.nii.gz", "central.nii.gz", "outer.nii.gz"):
                _save(output_dir / name, np.ones(shape, dtype=np.float32), affine)
            (output_dir / "provenance.json").write_text('{"fixture":true}\n')
        elif stage == "segmentation-refine":
            source = nib.load(command[command.index("--segmentation") + 1])
            _save(Path(command[command.index("--output") + 1]), np.asanyarray(source.dataobj).astype(np.int16), affine)
            Path(command[command.index("--summary") + 1]).write_text('{"source_faithful":true}\n')
        elif stage == "base-surface-thickness":
            mesh = SurfaceMesh(
                vertices=np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
                faces=np.asarray([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]], dtype=np.int64),
            )
            for name in ("inner_surface.vtk", "central_surface.vtk", "outer_surface.vtk", "inner_surface_thickness.vtk"):
                write_legacy_vtk(mesh, output_dir / name)
            _save(output_dir / "thickness_levelset.nii.gz", np.zeros(shape, dtype=np.float32), affine)
            (output_dir / "thickness_summary.json").write_text('{"fixture":true}\n')
            (output_dir / "provenance.json").write_text('{"fixture":true}\n')
        else:
            raise AssertionError(f"unexpected fake CLI command: {command}")
        return 0

    monkeypatch.setattr(resume_module, "cli_main", fake_cli)
    return calls


def test_tail_imports_hashed_upstream_and_resumes_verified_stages(tmp_path, monkeypatch):
    source, maseg = _source_run(tmp_path)
    calls = _install_fake_cli(monkeypatch, tmp_path / "source-t1.nii.gz")
    output = tmp_path / "tail-run"
    manifest_path = run_base_resume_from_toads(source_run=source, maseg_path=maseg, output_dir=output)
    manifest = json.loads(manifest_path.read_text())
    assert manifest["workflow_profile"] == "base-resume-from-toads"
    assert manifest["status"] == "complete"
    assert manifest["source_toads_second"]["source_stage"] == "toads-second"
    assert manifest["source_toads_second"]["source_stage_outputs"]
    assert set(manifest["stages"]) == {
        "base-tissue-preparation", "cruise-core", "maseg-refinement", "surface-thickness"
    }
    assert all(item["status"] == "complete" and item["dependency_signature"]
               for item in manifest["stages"].values())
    prep_command = manifest["stages"]["base-tissue-preparation"]["command"]
    assert prep_command[prep_command.index("--edge-width") + 1] == "5"
    core_command = manifest["stages"]["cruise-core"]["command"]
    assert core_command[core_command.index("--ace-diffusion-iterations") + 1] == "25"
    assert core_command[core_command.index("--tgdm-initial-iterations") + 1] == "6"
    assert core_command[core_command.index("--tgdm-inner-iterations") + 1] == "4"
    assert core_command[core_command.index("--tgdm-central-iterations") + 1] == "7"
    assert "--source-faithful" in manifest["stages"]["maseg-refinement"]["command"]
    count = len(calls)
    assert run_base_resume_from_toads(
        source_run=source, maseg_path=maseg, output_dir=output, resume=True
    ) == manifest_path
    assert len(calls) == count


def test_tail_rejects_changed_imported_toads_hash_before_creating_output(tmp_path, monkeypatch):
    source, maseg = _source_run(tmp_path)
    _install_fake_cli(monkeypatch, tmp_path / "source-t1.nii.gz")
    source_manifest = json.loads(source.read_text())
    gm_path = Path(next(item["path"] for item in source_manifest["stages"]["toads-second"]["outputs"]
                        if Path(item["path"]).name == "stripped_t1_gm.nii.gz"))
    with gm_path.open("ab") as stream:
        stream.write(b"tamper")
    output = tmp_path / "must-not-exist"
    with pytest.raises(ValueError, match="hash/size changed"):
        run_base_resume_from_toads(source_run=source, maseg_path=maseg, output_dir=output)
    assert not output.exists()


def test_implementation_change_invalidates_tail_and_reruns_downstream(tmp_path, monkeypatch):
    source, maseg = _source_run(tmp_path)
    calls = _install_fake_cli(monkeypatch, tmp_path / "source-t1.nii.gz")
    output = tmp_path / "tail-run"
    run_base_resume_from_toads(source_run=source, maseg_path=maseg, output_dir=output)
    assert len(calls) == 4
    monkeypatch.setattr(resume_module, "_implementation_signature", lambda: "d" * 64)
    run_base_resume_from_toads(source_run=source, maseg_path=maseg, output_dir=output, resume=True)
    assert len(calls) == 8
    manifest = json.loads((output / "run.json").read_text())
    assert manifest["implementation_signature"] == "d" * 64
    assert all(manifest["stages"][name]["attempt"] == 2 for name in manifest["stages"])


def test_tail_rejects_wrong_grid_maseg(tmp_path):
    source, maseg = _source_run(tmp_path)
    wrong = tmp_path / "wrong-grid.nii.gz"
    image = nib.load(str(maseg))
    affine = image.affine.copy()
    affine[0, 3] += 3.0
    _save(wrong, np.asanyarray(image.dataobj), affine)
    with pytest.raises(ValueError, match="must match the source run T1 grid"):
        run_base_resume_from_toads(source_run=source, maseg_path=wrong, output_dir=tmp_path / "out")
