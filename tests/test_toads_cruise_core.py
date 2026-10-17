from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.cli import main
from pymacruise.cropping import crop_image
from pymacruise.tgdm import TGDMBaselineParameters
import pymacruise.toads_cruise_core as core_module
from pymacruise.toads_cruise_core import run_toads_cruise_core


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "toads_tissue_crop"


def _inputs() -> tuple[nib.spatialimages.SpatialImage, ...]:
    return tuple(
        nib.load(str(FIXTURE / filename))
        for filename in ("gm.nii.gz", "wmfill.nii.gz", "wmmask.nii.gz")
    )


def _run_core():
    return run_toads_cruise_core(
        *_inputs(),
        crop_edge_width=2,
        ace_diffusion_iterations=1,
        gvf_iterations=200,
        tgdm_parameters=TGDMBaselineParameters(
            white_matter_iso=0.6,
            initial_iterations=1,
            inner_iterations=1,
            central_shrink_iterations=1,
            central_iterations=1,
            outer_iterations=1,
        ),
    )


def test_cruise_core_runs_ace_gvf_tgdm_and_restores_native_ras_grid() -> None:
    result = _run_core()
    reference = nib.load(str(FIXTURE / "gm.nii.gz"))
    assert result.preparation.stage_order == (
        "cropByPositiveGrayMatter",
        "topologyCorrectWhiteMatterFill",
    )
    assert result.ace.enhanced_gray_matter.shape == (11, 11, 9)
    assert result.gvf.vector_field.shape == (11, 11, 9, 3)
    assert result.tgdm_crop.inner.shape == (11, 11, 9)
    for image in (result.inner_native, result.central_native, result.outer_native):
        assert image.shape == reference.shape
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(image.affine, reference.affine)
        assert image.get_data_dtype() == np.dtype(np.float32)
    assert np.all(np.asarray(result.inner_native.dataobj)[0] == 20.0)


def test_cruise_core_crops_and_preserves_brain_mask_level_set_for_tgdm(monkeypatch) -> None:
    original_tgdm = core_module.tgdm_nested_images
    observed: dict[str, nib.spatialimages.SpatialImage | None] = {}

    def capture_tgdm(*args, **kwargs):
        observed["brain_mask_level_set"] = kwargs.get("brain_mask_level_set")
        return original_tgdm(*args, **kwargs)

    monkeypatch.setattr(core_module, "tgdm_nested_images", capture_tgdm)
    gm, wm, mask = _inputs()
    x, y, z = np.indices(gm.shape, dtype=np.float32)
    level_set = ((x - 6.0) + 0.25 * (y - 5.0) - 0.1 * z).astype(np.float32)
    source_level_set = nib.Nifti1Image(level_set, gm.affine)
    result = run_toads_cruise_core(
        gm, wm, mask, brain_mask_level_set=source_level_set,
        crop_edge_width=2, ace_diffusion_iterations=1, gvf_iterations=200,
        tgdm_parameters=TGDMBaselineParameters(
            initial_iterations=0, inner_iterations=0, central_shrink_iterations=0,
            central_iterations=0, outer_iterations=0,
        ),
    )

    observed_level_set = observed["brain_mask_level_set"]
    assert observed_level_set is not None
    expected = crop_image(source_level_set, result.preparation.crop_parameters)
    np.testing.assert_array_equal(observed_level_set.dataobj, expected.dataobj)
    np.testing.assert_allclose(observed_level_set.affine, result.preparation.gray_matter.affine)
    assert float(np.max(np.abs(np.asarray(observed_level_set.dataobj)))) > 1.0


def test_cruise_core_cli_writes_typed_native_grid_products_and_provenance(tmp_path) -> None:
    output_dir = tmp_path / "core"
    provenance = tmp_path / "core.json"
    assert main([
        "toads-cruise-core",
        "--gm", str(FIXTURE / "gm.nii.gz"),
        "--wm-fill", str(FIXTURE / "wmfill.nii.gz"),
        "--wm-mask", str(FIXTURE / "wmmask.nii.gz"),
        "--brain-mask-level-set", str(FIXTURE / "gm.nii.gz"),
        "--edge-width", "2",
        "--ace-diffusion-iterations", "1",
        "--gvf-iterations", "200",
        "--tgdm-initial-iterations", "1",
        "--tgdm-inner-iterations", "1",
        "--tgdm-central-iterations", "1",
        "--wm-iso-level", "0.43",
        "--use-updated-csf-in-ace", "0",
        "--output-dir", str(output_dir),
        "--provenance", str(provenance),
    ]) == 0
    assert provenance.is_file()
    parameters = json.loads(provenance.read_text())["parameters"]
    assert parameters["tgdm"]["white_matter_iso"] == 0.43
    assert parameters["ace"]["updated_multi_atlas_csf_mode"] == 0
    assert parameters["brain_mask_level_set"] == str((FIXTURE / "gm.nii.gz").resolve())
    reference = nib.load(str(FIXTURE / "gm.nii.gz"))
    expected = (
        "inner.nii.gz", "central.nii.gz", "outer.nii.gz",
    )
    for filename in expected:
        image = nib.load(str(output_dir / filename))
        assert image.shape == reference.shape
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert image.get_data_dtype() == np.dtype(np.float32)
    gvf = nib.load(str(output_dir / "gvf.nii.gz"))
    assert gvf.shape == (11, 11, 9, 3)
    assert gvf.get_data_dtype() == np.dtype(np.float32)


def test_cruise_core_default_iteration_path_runs_on_compact_synthetic_maps() -> None:
    result = run_toads_cruise_core(*_inputs(), crop_edge_width=2)
    assert result.ace.enhanced_gray_matter.shape == (11, 11, 9)
    assert result.gvf.vector_field.shape == (11, 11, 9, 3)
    assert result.tgdm_crop.execution.inner.completed_outer_passes > 0
    assert result.tgdm_crop.execution.central.completed_outer_passes > 0
    assert result.tgdm_crop.execution.outer.completed_outer_passes > 0
    for image in (result.inner_native, result.central_native, result.outer_native):
        values = np.asarray(image.dataobj)
        assert np.isfinite(values).all()
        assert image.shape == (13, 12, 11)


def test_cruise_core_rejects_boundary_free_ace_mode_two() -> None:
    with pytest.raises(ValueError, match="ACE mode 2 requires.*boundary_csf"):
        run_toads_cruise_core(
            *_inputs(),
            ace_use_updated_csf_in_ace=2,
        )
