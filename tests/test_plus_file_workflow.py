from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.image_io import save_nifti_ras
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


def test_plus_file_workflow_runs_resumes_and_repairs_invalid_tail(tmp_path: Path) -> None:
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

    with pytest.raises(ValueError, match="inputs or settings differ"):
        run_plus_file_workflow(**inputs, crop_edge_width=4, resume=True)
