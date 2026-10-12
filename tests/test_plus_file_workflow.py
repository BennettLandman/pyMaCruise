from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.image_io import save_nifti_ras
from pymacruise.boundary import clean_csf_boundary, multi_atlas_csf_boundary
from pymacruise.cropping import (
    crop_image,
    crop_image_with_reference_bounds,
    find_crop_parameters,
)
from pymacruise.plus_file_workflow import run_plus_file_workflow
from pymacruise.surfaces import read_legacy_vtk


FIXTURE = Path(__file__).parent / "data/synthetic/plus_pre_tgdm_composition"


def _load(name: str) -> nib.Nifti1Image:
    image = nib.load(str(FIXTURE / name))
    assert isinstance(image, nib.Nifti1Image)
    return image


def _save_channel(tmp_path: Path, name: str, values: np.ndarray, reference: nib.Nifti1Image) -> Path:
    image = nib.Nifti1Image(np.asarray(values, dtype=np.float32), reference.affine, reference.header)
    return save_nifti_ras(image, tmp_path / f"{name}.nii.gz")


def test_plus_file_workflow_runs_resumes_and_repairs_invalid_tail(tmp_path: Path, capsys) -> None:
    segmentation = _load("segmentation.nii.gz")
    atlas_gm = _load("multi_atlas_gm.nii.gz")
    atlas_wm = _load("multi_atlas_wm.nii.gz")
    atlas_mask = _load("multi_atlas_brain_mask.nii.gz")
    toads_gm = np.asarray(_load("toads_gm.nii.gz").dataobj, dtype=np.float32)
    toads_wm = np.asarray(_load("toads_wm.nii.gz").dataobj, dtype=np.float32)
    brain = _load("multi_atlas_brain_mask.nii.gz")
    channels = {
        "csf": np.zeros(segmentation.shape, dtype=np.float32),
        "gm": toads_gm,
        "wm_fill": toads_wm,
        "wm_mask": np.ones(segmentation.shape, dtype=np.float32),
    }
    inputs = {
        "t1_reference_path": _save_channel(tmp_path, "t1", np.zeros(segmentation.shape), segmentation),
        "segmentation_path": FIXTURE / "segmentation.nii.gz",
        "dura_csf_path": _save_channel(tmp_path, "dura_csf", channels["csf"], segmentation),
        "dura_gm_path": _save_channel(tmp_path, "dura_gm", channels["gm"], segmentation),
        "dura_wm_fill_path": _save_channel(tmp_path, "dura_wm_fill", channels["wm_fill"], segmentation),
        "dura_wm_mask_path": _save_channel(tmp_path, "dura_wm_mask", channels["wm_mask"], segmentation),
        "multi_atlas_gm_path": FIXTURE / "multi_atlas_gm.nii.gz",
        "multi_atlas_wm_path": FIXTURE / "multi_atlas_wm.nii.gz",
        "multi_atlas_mask_path": FIXTURE / "multi_atlas_brain_mask.nii.gz",
        "brain_mask_path": FIXTURE / "multi_atlas_brain_mask.nii.gz",
        "output_dir": tmp_path / "plus-run",
    }

    mismatched_mask = save_nifti_ras(
        nib.Nifti1Image(
            np.zeros((segmentation.shape[0] - 1, *segmentation.shape[1:]), dtype=np.uint8),
            segmentation.affine,
        ),
        tmp_path / "mismatched-brain-mask.nii.gz",
    )
    with pytest.raises(ValueError, match="brain_mask must match the T1 reference grid"):
        run_plus_file_workflow(
            **{**inputs, "brain_mask_path": mismatched_mask},
            crop_edge_width=5,
        )

    first = run_plus_file_workflow(**inputs, crop_edge_width=5)
    assert first.resumed_stages == ()
    manifest = json.loads(first.manifest_path.read_text())
    assert manifest["status"] == "complete"
    assert set(manifest["stages"]) == {"transfer_crop", "topology_ace_gvf_tgdm", "surfaces_thickness"}
    assert set(first.outputs) == {
        "inner_level_set", "central_level_set", "outer_level_set", "thickness_volume",
        "inner_surface", "central_surface", "outer_surface", "thickness_metrics",
    }
    for name, path in first.outputs.items():
        assert path.is_file()
        if path.suffix == ".gz":
            image = nib.load(str(path))
            assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
            assert image.shape == segmentation.shape
            assert np.isfinite(np.asarray(image.dataobj)).all()
        elif path.suffix == ".vtk":
            mesh = read_legacy_vtk(path)
            assert mesh.coordinate_frame == "RAS"
            assert mesh.units == "mm"
            assert np.isfinite(mesh.vertices).all()
    metrics = json.loads(first.outputs["thickness_metrics"].read_text())
    assert metrics["spatial_convention"] == "RAS"
    assert metrics["surface_units"] == "mm"

    # Exercise the public CLI against the same files, resuming the valid run.
    from pymacruise.cli import main

    cli_flags = {
        "t1_reference_path": "--t1-reference",
        "segmentation_path": "--segmentation",
        "dura_csf_path": "--dura-csf",
        "dura_gm_path": "--dura-gm",
        "dura_wm_fill_path": "--dura-wm-fill",
        "dura_wm_mask_path": "--dura-wm-mask",
        "multi_atlas_gm_path": "--multi-atlas-gm",
        "multi_atlas_wm_path": "--multi-atlas-wm",
        "multi_atlas_mask_path": "--multi-atlas-mask",
        "brain_mask_path": "--brain-mask",
        "output_dir": "--output-dir",
    }
    cli_arguments = ["plus-run"]
    for name, flag in cli_flags.items():
        cli_arguments.extend((flag, str(inputs[name])))
    cli_arguments.extend(("--resume", "--crop-edge-width", "5"))
    assert main(cli_arguments) == 0
    assert "Verified stages skipped: transfer_crop, topology_ace_gvf_tgdm, surfaces_thickness" in capsys.readouterr().out

    resumed = run_plus_file_workflow(**inputs, crop_edge_width=5, resume=True)
    assert resumed.resumed_stages == ("transfer_crop", "topology_ace_gvf_tgdm", "surfaces_thickness")

    resumed.outputs["thickness_volume"].write_bytes(b"corrupt")
    repaired = run_plus_file_workflow(**inputs, crop_edge_width=5, resume=True)
    assert repaired.resumed_stages == ("transfer_crop", "topology_ace_gvf_tgdm")
    assert np.isfinite(np.asarray(nib.load(str(repaired.outputs["thickness_volume"])).dataobj)).all()

    manifest = json.loads(repaired.manifest_path.read_text())
    manifest["status"] = "complete"
    manifest["stages"]["surfaces_thickness"]["status"] = "running"
    repaired.manifest_path.write_text(json.dumps(manifest))
    interrupted = run_plus_file_workflow(**inputs, crop_edge_width=5, resume=True)
    assert interrupted.resumed_stages == ("transfer_crop", "topology_ace_gvf_tgdm")
    assert json.loads(interrupted.manifest_path.read_text())["status"] == "complete"

    changed_settings = run_plus_file_workflow(**inputs, crop_edge_width=4, resume=True)
    assert changed_settings.resumed_stages == ()
    changed_manifest = json.loads(changed_settings.manifest_path.read_text())
    assert changed_manifest["status"] == "complete"
    assert changed_manifest["previous_run_signature"] != changed_manifest["run_signature"]
    assert changed_manifest["settings"]["crop_edge_width"] == 4


def test_plus_file_workflow_jist_unchanged_mode_crops_smaller_dura_extents(tmp_path: Path) -> None:
    segmentation = _load("segmentation.nii.gz")
    shape = tuple(int(size - 2) for size in segmentation.shape)
    source_gm = np.asarray(_load("toads_gm.nii.gz").dataobj, dtype=np.float32)[: shape[0], : shape[1], : shape[2]]
    source_wm = np.asarray(_load("toads_wm.nii.gz").dataobj, dtype=np.float32)[: shape[0], : shape[1], : shape[2]]
    source_mask = np.ones(shape, dtype=np.float32)
    source_csf = np.zeros(shape, dtype=np.float32)
    inputs = {
        "t1_reference_path": _save_channel(tmp_path, "jist-t1", np.zeros(segmentation.shape), segmentation),
        "segmentation_path": FIXTURE / "segmentation.nii.gz",
        "dura_csf_path": _save_channel(tmp_path, "jist-dura-csf", source_csf, segmentation),
        "dura_gm_path": _save_channel(tmp_path, "jist-dura-gm", source_gm, segmentation),
        "dura_wm_fill_path": _save_channel(tmp_path, "jist-dura-wm-fill", source_wm, segmentation),
        "dura_wm_mask_path": _save_channel(tmp_path, "jist-dura-wm-mask", source_mask, segmentation),
        "multi_atlas_gm_path": FIXTURE / "multi_atlas_gm.nii.gz",
        "multi_atlas_wm_path": FIXTURE / "multi_atlas_wm.nii.gz",
        "multi_atlas_mask_path": FIXTURE / "multi_atlas_brain_mask.nii.gz",
        "brain_mask_path": FIXTURE / "multi_atlas_brain_mask.nii.gz",
        "output_dir": tmp_path / "jist-unchanged-run",
        "dura_transfer_mode": "jist-unchanged",
    }

    run = run_plus_file_workflow(**inputs, crop_edge_width=1)
    manifest = json.loads(run.manifest_path.read_text())
    assert manifest["settings"]["dura_transfer_mode"] == "jist-unchanged"
    crop_attempt = Path(manifest["stages"]["transfer_crop"]["outputs"]["toads_gm"]["path"]).parent
    crop_parameters = find_crop_parameters(_load("multi_atlas_gm.nii.gz"), edge_width=1)
    expected_gm = crop_image_with_reference_bounds(
        nib.Nifti1Image(source_gm, segmentation.affine), crop_parameters
    )
    actual_gm = nib.load(str(crop_attempt / "toads_gm.nii.gz"))
    assert actual_gm.shape == expected_gm.shape
    np.testing.assert_array_equal(actual_gm.dataobj, expected_gm.dataobj)
    assert nib.aff2axcodes(actual_gm.affine) == ("R", "A", "S")
    assert set(run.outputs) == {
        "inner_level_set", "central_level_set", "outer_level_set", "thickness_volume",
        "inner_surface", "central_surface", "outer_surface", "thickness_metrics",
    }


def test_plus_file_workflow_uses_toads_wm_mask_for_boundary_cleanup(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The selected source cleans the MaSeg boundary with TOADS WM mask (>0.01)."""
    import pymacruise.plus_file_workflow as workflow

    segmentation = _load("segmentation.nii.gz")
    atlas_gm = _load("multi_atlas_gm.nii.gz")
    crop = find_crop_parameters(atlas_gm, threshold=0.0, edge_width=0)
    cropped_segmentation = crop_image(segmentation, crop)
    raw_boundary = multi_atlas_csf_boundary(np.asarray(cropped_segmentation.dataobj))
    boundary_voxels = np.argwhere(raw_boundary > 0)
    assert boundary_voxels.size > 0, "fixture must exercise nonempty source boundary cleanup"

    # Equality to 0.01 must survive the source's strict `> 0.01` test. Put a
    # slightly larger mask value on one real boundary voxel so both branches
    # of the threshold are observable in one workflow execution.
    wm_mask = np.full(segmentation.shape, 0.01, dtype=np.float32)
    full_boundary_voxel = boundary_voxels[0] + np.asarray(crop.offset)
    wm_mask[tuple(full_boundary_voxel)] = np.float32(0.0101)
    wm_fill = np.full(segmentation.shape, 255.0, dtype=np.float32)

    toads_gm = np.asarray(_load("toads_gm.nii.gz").dataobj, dtype=np.float32)
    channels = {
        "csf": np.zeros(segmentation.shape, dtype=np.float32),
        "gm": toads_gm,
        "wm_fill": wm_fill,
        "wm_mask": wm_mask,
    }
    inputs = {
        "t1_reference_path": _save_channel(tmp_path, "boundary-t1", np.zeros(segmentation.shape), segmentation),
        "segmentation_path": FIXTURE / "segmentation.nii.gz",
        "dura_csf_path": _save_channel(tmp_path, "boundary-dura-csf", channels["csf"], segmentation),
        "dura_gm_path": _save_channel(tmp_path, "boundary-dura-gm", channels["gm"], segmentation),
        "dura_wm_fill_path": _save_channel(tmp_path, "boundary-dura-wm-fill", channels["wm_fill"], segmentation),
        "dura_wm_mask_path": _save_channel(tmp_path, "boundary-dura-wm-mask", channels["wm_mask"], segmentation),
        "multi_atlas_gm_path": FIXTURE / "multi_atlas_gm.nii.gz",
        "multi_atlas_wm_path": FIXTURE / "multi_atlas_wm.nii.gz",
        "multi_atlas_mask_path": FIXTURE / "multi_atlas_brain_mask.nii.gz",
        "brain_mask_path": FIXTURE / "multi_atlas_brain_mask.nii.gz",
        "output_dir": tmp_path / "boundary-run",
    }

    observed: dict[str, nib.spatialimages.SpatialImage] = {}
    original = workflow.run_plus_topology_and_tgdm

    def traced_topology(seg, wm_for_boundary, *args, **kwargs):
        observed["segmentation"] = seg
        observed["wm_for_boundary"] = wm_for_boundary
        observed["result"] = original(seg, wm_for_boundary, *args, **kwargs)
        return observed["result"]

    monkeypatch.setattr(workflow, "run_plus_topology_and_tgdm", traced_topology)
    run_plus_file_workflow(**inputs, crop_edge_width=0)

    result = observed["result"]
    assert int(np.asarray(result.boundary.raw_boundary.dataobj).sum()) > 0
    actual_mask = np.asarray(observed["wm_for_boundary"].dataobj)
    expected_cleaned = clean_csf_boundary(
        np.asarray(result.boundary.raw_boundary.dataobj), actual_mask, mode="toads_wm"
    )
    np.testing.assert_array_equal(result.boundary.cleaned_boundary.dataobj, expected_cleaned)

    local_equal_voxels = (np.asarray(result.boundary.raw_boundary.dataobj) > 0) & np.isclose(
        actual_mask, 0.01, rtol=0, atol=1e-8
    )
    local_above_voxels = (np.asarray(result.boundary.raw_boundary.dataobj) > 0) & (actual_mask > 0.01)
    assert np.any(local_equal_voxels), "boundary voxels at exactly 0.01 must remain"
    assert np.any(local_above_voxels), "fixture must include a boundary voxel above 0.01"
    assert np.all(np.asarray(result.boundary.cleaned_boundary.dataobj)[local_equal_voxels] > 0)
    assert np.all(np.asarray(result.boundary.cleaned_boundary.dataobj)[local_above_voxels] == 0)

    # WM fill is intentionally different from the routed mask. Confirm the
    # resulting cleanup differs too, so the regression catches routing fill.
    cropped_fill_path = json.loads((inputs["output_dir"] / "run.json").read_text())["stages"]["transfer_crop"]["outputs"]["toads_wm_fill"]["path"]
    cropped_fill = np.asarray(nib.load(cropped_fill_path).dataobj)
    assert not np.array_equal(cropped_fill, actual_mask)
    wrong_fill_result = clean_csf_boundary(
        np.asarray(result.boundary.raw_boundary.dataobj), cropped_fill, mode="toads_wm"
    )
    assert not np.array_equal(wrong_fill_result, result.boundary.cleaned_boundary.dataobj)
