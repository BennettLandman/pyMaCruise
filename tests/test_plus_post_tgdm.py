from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

import pymacruise.plus_pipeline as plus_pipeline
from pymacruise.cropping import CropParameters
from pymacruise.plus_pipeline import (
    extract_plus_surfaces,
    run_plus_post_tgdm_stages,
    run_plus_thickness_stage,
    uncrop_plus_tgdm_outputs,
)
from pymacruise.tgdm import TGDMMaCRUISEImageResult


FIXTURE = Path(__file__).parent / "data/synthetic/plus_post_tgdm"
FIELDS = (
    "initial_level_set",
    "tissue_sign",
    "initial_inner",
    "inner",
    "central",
    "outer",
)


def _cropped_tgdm() -> TGDMMaCRUISEImageResult:
    images = {name: nib.load(FIXTURE / f"{name}.nii.gz") for name in FIELDS}
    # Execution provenance is intentionally absent: these are fixed synthetic
    # upstream stage products, not outputs of the full TGDM evolution.
    return TGDMMaCRUISEImageResult(
        **images,
        execution=None,  # type: ignore[arg-type]
    )


def _crop_parameters() -> CropParameters:
    reference = nib.load(FIXTURE / "outer.nii.gz")
    affine = reference.affine.copy()
    translation = np.eye(4)
    translation[:3, 3] = (-2, -2, -2)
    original_affine = affine @ translation
    return CropParameters(
        original_shape=(17, 17, 17),
        bounds=((2, 14), (2, 14), (2, 14)),
        background=0.0,
        original_affine=original_affine,
    )


def test_plus_post_tgdm_composes_public_stages_in_source_order(monkeypatch) -> None:
    events: list[str] = []
    direct_uncrop = plus_pipeline.uncrop_plus_tgdm_outputs
    direct_surfaces = plus_pipeline.extract_plus_surfaces
    direct_thickness = plus_pipeline.run_plus_thickness_stage

    def traced_uncrop(*args, **kwargs):
        events.append("uncrop")
        return direct_uncrop(*args, **kwargs)

    def traced_surfaces(*args, **kwargs):
        events.append("surfaces")
        return direct_surfaces(*args, **kwargs)

    def traced_thickness(*args, **kwargs):
        events.append("thickness")
        return direct_thickness(*args, **kwargs)

    monkeypatch.setattr(plus_pipeline, "uncrop_plus_tgdm_outputs", traced_uncrop)
    monkeypatch.setattr(plus_pipeline, "extract_plus_surfaces", traced_surfaces)
    monkeypatch.setattr(plus_pipeline, "run_plus_thickness_stage", traced_thickness)

    cropped = _cropped_tgdm()
    parameters = _crop_parameters()
    saved_inputs = {
        name: (
            np.asarray(getattr(cropped, name).dataobj).copy(),
            getattr(cropped, name).affine.copy(),
        )
        for name in FIELDS
    }
    actual = run_plus_post_tgdm_stages(cropped, parameters)

    assert events == ["uncrop", "surfaces", "thickness"]

    expected_level_sets = uncrop_plus_tgdm_outputs(cropped, parameters)
    expected_surfaces = extract_plus_surfaces(expected_level_sets)
    expected_thickness = run_plus_thickness_stage(
        expected_level_sets, expected_surfaces
    )
    for name in FIELDS:
        observed = getattr(actual.level_sets, name)
        expected = getattr(expected_level_sets, name)
        np.testing.assert_array_equal(observed.dataobj, expected.dataobj)
        np.testing.assert_array_equal(observed.affine, expected.affine)
        assert observed.shape == (17, 17, 17)
        assert nib.aff2axcodes(observed.affine) == ("R", "A", "S")

    for name in ("inner", "central", "outer"):
        observed = getattr(actual.surfaces, name)
        expected = getattr(expected_surfaces, name)
        np.testing.assert_array_equal(observed.vertices, expected.vertices)
        np.testing.assert_array_equal(observed.faces, expected.faces)
        assert observed.coordinate_frame == "RAS"
        assert observed.units == "mm"

    np.testing.assert_array_equal(
        actual.thickness.volume.image.dataobj,
        expected_thickness.volume.image.dataobj,
    )
    assert actual.thickness.volume.mean_mm == expected_thickness.volume.mean_mm
    np.testing.assert_array_equal(
        actual.thickness.surface.inner_surface.point_data["mean_thickness_mm"],
        expected_thickness.surface.inner_surface.point_data["mean_thickness_mm"],
    )
    assert (
        actual.thickness.surface.mean_mm
        == expected_thickness.surface.mean_mm
    )
    for name, (saved_data, saved_affine) in saved_inputs.items():
        observed = getattr(cropped, name)
        np.testing.assert_array_equal(observed.dataobj, saved_data)
        np.testing.assert_array_equal(observed.affine, saved_affine)


def test_plus_post_tgdm_keeps_selected_uncrop_boundary_and_mesh_risk_visible() -> None:
    actual = run_plus_post_tgdm_stages(_cropped_tgdm(), _crop_parameters())
    # JIST's selected exclusive maximum leaves the final source plane at the
    # 20-valued MaTGDMcrop exterior. This fixture pins that behavior at the
    # crop edge; zero would be changed to -0.1 by the selected surface adapter.
    assert np.all(np.asarray(actual.level_sets.outer.dataobj)[14, :, :] == 20.0)
    assert all(
        nib.aff2axcodes(getattr(actual.level_sets, name).affine) == ("R", "A", "S")
        for name in FIELDS
    )
