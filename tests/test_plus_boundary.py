from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.boundary import clean_csf_boundary_image, multi_atlas_csf_boundary_image
from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_boundary_stage


FIXTURE = Path(__file__).parent / "data/synthetic/plus_boundary"


def test_plus_boundary_stage_composes_configured_3d_boundary_and_toads_wm_cleanup() -> None:
    segmentation = nib.load(FIXTURE / "segmentation.nii.gz")
    toads_wm = nib.load(FIXTURE / "toads_wm_mask.nii.gz")
    options = PlusPipelineOptions(boundary_overlook_labels=(120,))
    actual = run_plus_boundary_stage(segmentation, toads_wm, options)
    raw = multi_atlas_csf_boundary_image(
        segmentation, overlooked_labels=options.boundary_overlook_labels
    )
    cleaned = clean_csf_boundary_image(raw, toads_wm, mode="toads_wm")

    np.testing.assert_array_equal(actual.raw_boundary.dataobj, raw.dataobj)
    np.testing.assert_array_equal(actual.cleaned_boundary.dataobj, cleaned.dataobj)
    assert nib.aff2axcodes(actual.raw_boundary.affine) == ("R", "A", "S")
    assert np.asarray(actual.raw_boundary.dataobj).any()
    assert np.any(
        np.asarray(actual.raw_boundary.dataobj)
        != np.asarray(actual.cleaned_boundary.dataobj)
    )
    boundary_without_overlook = multi_atlas_csf_boundary_image(
        segmentation, overlooked_labels=()
    )
    assert np.any(
        np.asarray(actual.raw_boundary.dataobj)
        != np.asarray(boundary_without_overlook.dataobj)
    )
    wm_voxels = np.asarray(toads_wm.dataobj) > 0.01
    assert np.all(np.asarray(actual.cleaned_boundary.dataobj)[wm_voxels] == 0)


def test_plus_boundary_stage_rejects_nonselected_erosion_method() -> None:
    segmentation = nib.load(FIXTURE / "segmentation.nii.gz")
    toads_wm = nib.load(FIXTURE / "toads_wm_mask.nii.gz")
    options = PlusPipelineOptions(boundary_erosion_method="3Dmatlab")
    with pytest.raises(ValueError, match="only ErodeMethod='3D'"):
        run_plus_boundary_stage(segmentation, toads_wm, options)
