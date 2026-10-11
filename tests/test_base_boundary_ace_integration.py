"""Focused checks for the selected base MaCRUISE boundary-to-ACE route.

The source V3 route builds a MaSeg boundary on the shared crop, removes its
TOADS-WM support and then its updated-WM support, and supplies that map to ACE
with ``UseUpdatedCSFinACE=2``. These checks intentionally describe the data
contract the restartable workflow must wire.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np

from pymacruise import base_tissue_preparation as preparation_module
from pymacruise import cli as cli_module
from pymacruise.base_tissue_preparation import prepare_base_tissue_maps
from pymacruise.boundary import (
    clean_csf_boundary_image,
    multi_atlas_csf_boundary,
)
from pymacruise.cli import main
from pymacruise.tgdm import TGDMBaselineParameters
from pymacruise.toads_cruise_core import run_toads_cruise_core


def _preparation_inputs() -> tuple[nib.Nifti1Image, ...]:
    """Small source-label phantom with boundary and both WM-cleanup cases."""
    shape = (11, 11, 11)
    affine = np.eye(4)
    labels = np.zeros(shape, dtype=np.int16)
    labels[2:9, 2:9, 2:9] = 100
    labels[4:7, 4:7, 4:7] = 4
    labels[3:5, 4:7, 4:7] = 101
    labels[7, 5, 5] = 120  # ignored by source MinFilter3D
    gm = np.full(shape, 0.5, dtype=np.float32)
    wm = np.full(shape, 0.0, dtype=np.float32)
    toads_mask = np.zeros(shape, dtype=np.uint8)
    # A TOADS WM voxel removes an otherwise-positive boundary point.
    toads_mask[5, 4, 5] = 1
    # Updated WM uses the source byte-domain threshold (>25.5).
    wm[5, 5, 5] = 255.0
    return tuple(
        nib.Nifti1Image(data, affine)
        for data in (labels, gm, wm, toads_mask)
    )


def test_base_preparation_emits_twice_cleaned_boundary_on_shared_crop() -> None:
    maseg, toads_gm, toads_wm, toads_wm_mask = _preparation_inputs()
    prepared = prepare_base_tissue_maps(maseg, toads_gm, toads_wm, toads_wm_mask)

    # Expected sequence mirrors run_GetMABoundary -> CleanBoundaryByToads ->
    # CleanBoundaryByUpdateWM from the selected V3 base MATLAB wrapper.
    cropped_labels = preparation_module.crop_image(
        maseg, prepared.crop_parameters
    )
    raw = multi_atlas_csf_boundary(np.asarray(cropped_labels.dataobj))
    raw_image = nib.Nifti1Image(raw, cropped_labels.affine)
    after_toads = clean_csf_boundary_image(
        raw_image, prepared.cropped_toads_mask, mode="toads_wm"
    )
    expected = clean_csf_boundary_image(
        after_toads, prepared.white_matter, mode="updated_wm"
    )

    actual_image = prepared.boundary_csf
    np.testing.assert_array_equal(
        np.asarray(actual_image.dataobj), np.asarray(expected.dataobj)
    )
    assert actual_image.shape == prepared.gray_matter.shape
    np.testing.assert_allclose(actual_image.affine, prepared.gray_matter.affine)
    assert nib.aff2axcodes(actual_image.affine) == ("R", "A", "S")
    assert prepared.provenance["boundary_csf"]["cleanup_order"] == [
        "toads_wm_gt_0.01",
        "updated_wm_gt_25.5",
    ]


def test_prepared_boundary_is_passed_into_ace_at_mode_two(monkeypatch) -> None:
    """Exercise preparation→restartable core→ACE with the same cropped map."""
    maseg, toads_gm, toads_wm, toads_wm_mask = _preparation_inputs()
    prepared = prepare_base_tissue_maps(maseg, toads_gm, toads_wm, toads_wm_mask)
    from pymacruise import toads_cruise_core as core_module

    real_ace = core_module.ace_enhance
    captured: dict[str, object] = {}

    def spy_ace(gray, white, **kwargs):
        captured.update(kwargs)
        return real_ace(gray, white, **kwargs)

    monkeypatch.setattr(core_module, "ace_enhance", spy_ace)
    run_toads_cruise_core(
        prepared.gray_matter,
        prepared.white_matter,
        prepared.mask,
        boundary_csf=prepared.boundary_csf,
        ace_use_updated_csf_in_ace=2,
        ace_diffusion_iterations=1,
        gvf_iterations=2,
        tgdm_parameters=TGDMBaselineParameters(
            white_matter_iso=0.5,
            initial_iterations=1,
            inner_iterations=1,
            central_shrink_iterations=1,
            central_iterations=1,
            outer_iterations=1,
        ),
        prepared_crop_parameters=prepared.crop_parameters,
    )

    passed_boundary = captured["multi_atlas_csf"]
    np.testing.assert_array_equal(
        np.asarray(passed_boundary.dataobj), np.asarray(prepared.boundary_csf.dataobj)
    )
    assert captured["use_updated_csf_in_ace"] == 2


def test_core_cli_routes_boundary_to_core_and_provenance(tmp_path: Path, monkeypatch) -> None:
    """The public core CLI passes its boundary file and records its provenance."""
    fixture = Path(__file__).parent / "data" / "synthetic" / "toads_tissue_crop"
    boundary = nib.load(str(fixture / "wmmask.nii.gz"))
    boundary_path = tmp_path / "boundary-csf.nii.gz"
    nib.save(boundary, boundary_path)
    seen: dict[str, object] = {}

    def fake_core(gm, wm_fill, wm_mask, **kwargs):
        seen.update(kwargs)
        image = nib.load(str(fixture / "gm.nii.gz"))
        zero = nib.Nifti1Image(np.zeros(image.shape, dtype=np.float32), image.affine)
        vectors = nib.Nifti1Image(np.zeros((*image.shape, 3), dtype=np.float32), image.affine)
        return SimpleNamespace(
            preparation=SimpleNamespace(
                gray_matter=image, white_matter_fill=image, white_matter_mask=image,
                topology_corrected_white_matter_fill=image,
                crop_parameters=SimpleNamespace(bounds=((0, 1),) * 3),
                removed_voxel_counts=None, stage_order=("boundary-preparation",),
            ),
            ace=SimpleNamespace(enhanced_gray_matter=image, skeleton=zero, thinned_skeleton=zero),
            gvf=SimpleNamespace(vector_field=vectors),
            inner_native=zero, central_native=zero, outer_native=zero,
        )

    monkeypatch.setattr(cli_module, "run_toads_cruise_core", fake_core)
    output = tmp_path / "core"
    provenance_path = tmp_path / "core-provenance.json"
    assert main([
        "toads-cruise-core", "--gm", str(fixture / "gm.nii.gz"),
        "--wm-fill", str(fixture / "wmfill.nii.gz"),
        "--wm-mask", str(fixture / "wmmask.nii.gz"),
        "--boundary-csf", str(boundary_path),
        "--use-updated-csf-in-ace", "2", "--output-dir", str(output),
        "--provenance", str(provenance_path),
    ]) == 0

    assert isinstance(seen.get("boundary_csf"), nib.spatialimages.SpatialImage)
    np.testing.assert_array_equal(seen["boundary_csf"].dataobj, boundary.dataobj)
    assert seen["ace_use_updated_csf_in_ace"] == 2
    record = json.loads(provenance_path.read_text())
    assert record["inputs"]["boundary_csf"] == str(boundary_path.resolve())
