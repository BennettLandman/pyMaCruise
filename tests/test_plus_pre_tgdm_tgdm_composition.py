from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

import pymacruise.plus_composition as composition
from pymacruise.plus_composition import run_plus_pre_tgdm_and_tgdm
from pymacruise.plus_pipeline import PlusPipelineOptions


PRE_FIXTURE = Path(__file__).parent / "data/synthetic/plus_pre_tgdm_composition"
TGDM_FIXTURE = Path(__file__).parent / "data/synthetic/plus_tgdm_composition"
PRE_NAMES = (
    "segmentation.nii.gz",
    "toads_wm_for_boundary.nii.gz",
    "multi_atlas_gm.nii.gz",
    "multi_atlas_wm.nii.gz",
    "multi_atlas_brain_mask.nii.gz",
    "toads_gm.nii.gz",
    "toads_wm.nii.gz",
    "toads_mask.nii.gz",
)
SEED_NAME = "topology_corrected_initial_membership.nii.gz"


def _inputs() -> tuple[nib.Nifti1Image, ...]:
    seed_provenance = json.loads((TGDM_FIXTURE / "provenance.json").read_text())
    seed_path = TGDM_FIXTURE / SEED_NAME
    assert hashlib.sha256(seed_path.read_bytes()).hexdigest() == seed_provenance[
        "sha256"
    ][SEED_NAME]
    inputs = tuple(nib.load(str(PRE_FIXTURE / name)) for name in PRE_NAMES) + (
        nib.load(str(seed_path)),
    )
    for image in inputs:
        assert image.shape == (25, 25, 25)
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        np.testing.assert_allclose(
            image.affine,
            np.diag((0.8, 0.8, 1.2, 1.0)),
            rtol=0.0,
            atol=1e-6,
        )
    return inputs


def test_plus_pre_tgdm_and_tgdm_routes_explicit_seed_and_preserves_grid(
    monkeypatch,
) -> None:
    events: list[str] = []
    routed: list[tuple[object, ...]] = []
    normalized_pre_inputs: list[tuple[object, ...]] = []
    original_pre = composition.run_plus_pre_tgdm_stages
    original_tgdm = composition.run_plus_tgdm_stage

    def traced_pre(*args, **kwargs):
        events.append("pre_tgdm")
        normalized_pre_inputs.append(args)
        return original_pre(*args, **kwargs)

    def traced_tgdm(*args, **kwargs):
        events.append("tgdm")
        routed.append(args)
        return original_tgdm(*args, **kwargs)

    monkeypatch.setattr(composition, "run_plus_pre_tgdm_stages", traced_pre)
    monkeypatch.setattr(composition, "run_plus_tgdm_stage", traced_tgdm)
    ras_inputs = _inputs()
    flipped_data = np.flip(np.asarray(ras_inputs[0].dataobj), axis=0).copy()
    left_affine = ras_inputs[0].affine.copy()
    left_affine[0, 0] *= -1
    left_affine[0, 3] = ras_inputs[0].affine[0, 0] * (ras_inputs[0].shape[0] - 1)
    non_ras_segmentation = nib.Nifti1Image(flipped_data, left_affine)
    assert nib.aff2axcodes(non_ras_segmentation.affine) == ("L", "A", "S")
    inputs = (non_ras_segmentation, *ras_inputs[1:])
    before = [np.asarray(image.dataobj).copy() for image in inputs]
    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        result = run_plus_pre_tgdm_and_tgdm(*inputs, options=PlusPipelineOptions())

    assert events == ["pre_tgdm", "tgdm"]
    assert nib.aff2axcodes(normalized_pre_inputs[0][0].affine) == ("R", "A", "S")
    np.testing.assert_array_equal(
        normalized_pre_inputs[0][0].dataobj,
        ras_inputs[0].dataobj,
    )
    assert len(routed) == 1
    assert routed[0][0] is inputs[8]
    assert routed[0][1] is result.pre_tgdm.tissue.wm
    assert routed[0][2] is result.pre_tgdm.tissue.gm
    assert routed[0][3] is result.pre_tgdm.ace.enhanced_gray_matter
    assert routed[0][4] is result.pre_tgdm.ace.thinned_skeleton
    assert routed[0][5] is result.pre_tgdm.gvf.vector_field

    for name in (
        "initial_level_set",
        "tissue_sign",
        "initial_inner",
        "inner",
        "central",
        "outer",
    ):
        output = getattr(result.tgdm, name)
        assert output.shape == inputs[8].shape
        np.testing.assert_allclose(output.affine, inputs[8].affine, rtol=0, atol=1e-6)
        assert nib.aff2axcodes(output.affine) == ("R", "A", "S")
        assert np.isfinite(np.asarray(output.dataobj)).all()
    for image, original in zip(inputs, before, strict=True):
        np.testing.assert_array_equal(image.dataobj, original)
