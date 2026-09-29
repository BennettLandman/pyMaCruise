from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_ace_stage


ROOT = Path(__file__).parent / "data/synthetic"
INPUTS = ROOT / "ace"
REFERENCE = ROOT / "ace_after_gaussian/java_reference"


def test_plus_ace_stage_uses_profile_options_and_matches_java_outputs() -> None:
    gray = nib.load(INPUTS / "gray_matter.nii.gz")
    white = nib.load(INPUTS / "white_matter.nii.gz")
    boundary = nib.load(INPUTS / "multi_atlas_csf.nii.gz")
    result = run_plus_ace_stage(gray, white, boundary, PlusPipelineOptions())

    assert result.enhanced_gray_matter.shape == gray.shape
    assert nib.aff2axcodes(result.enhanced_gray_matter.affine) == ("R", "A", "S")
    for image, filename in (
        (result.enhanced_gray_matter, "enhanced_gm_java.nii.gz"),
        (result.skeleton, "csf_skeleton_java.nii.gz"),
        (result.thinned_skeleton, "thinned_skeleton_java.nii.gz"),
    ):
        expected = np.asarray(nib.load(REFERENCE / filename).dataobj)
        actual = np.asarray(image.dataobj)
        if filename == "enhanced_gm_java.nii.gz":
            error = np.abs(actual.astype(np.float32) - expected.astype(np.float32))
            assert float(error.mean()) <= 1e-6
            assert float(error.max()) <= 5e-6
        else:
            np.testing.assert_array_equal(actual, expected)


def test_plus_ace_stage_rejects_unsupported_post_ace_update() -> None:
    gray = nib.load(INPUTS / "gray_matter.nii.gz")
    white = nib.load(INPUTS / "white_matter.nii.gz")
    boundary = nib.load(INPUTS / "multi_atlas_csf.nii.gz")
    options = PlusPipelineOptions(update_enhance=1)
    with pytest.raises(ValueError, match="only update_enhance=0"):
        run_plus_ace_stage(gray, white, boundary, options)
