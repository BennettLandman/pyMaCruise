from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.ace import (
    ace_anisotropic_diffusion,
    ace_enhance_arrays,
    ace_gaussian_blur,
    ace_seed_fields,
    ace_update_speed,
    _signed_distance_fmm,
)


ROOT = Path(__file__).parent / "data/synthetic"
INPUTS = ROOT / "ace"
REFERENCE = ROOT / "ace_after_gaussian/java_reference"


def _array(path: Path) -> np.ndarray:
    return np.asarray(nib.load(path).dataobj)


def _assert_close(actual: np.ndarray, expected: np.ndarray) -> None:
    error = np.abs(np.asarray(actual, dtype=np.float32) - np.asarray(expected, dtype=np.float32))
    assert float(error.mean()) <= 1.0e-6
    assert float(error.max()) <= 5.0e-6


def test_plus_after_gaussian_speed_and_ace_outputs_match_java() -> None:
    gray = _array(INPUTS / "gray_matter.nii.gz").astype(np.float32)
    white = _array(INPUTS / "white_matter.nii.gz").astype(np.float32)
    atlas = _array(INPUTS / "multi_atlas_csf.nii.gz")

    initial_phi = _signed_distance_fmm(np.float32(126.9) - white, 1000.0)
    seeds = ace_seed_fields(gray, white, initial_phi, intensity_threshold=126.9)
    guide = ace_gaussian_blur(initial_phi, 1.8)
    diffused = ace_anisotropic_diffusion(
        seeds.initial_csf_speed, guide, tensor_alpha=1.6, iterations=25
    )
    speed_before_merge = ace_gaussian_blur(1.0 - np.float32(0.9) * diffused, 1.8)
    expected_speed = _array(
        REFERENCE / "stages/11_5_initcsf_afterGaussian_mAcsf_java.nii.gz"
    )
    speed = ace_update_speed(
        diffused,
        multi_atlas_csf=atlas,
        gaussian_alpha=1.8,
        place_to_add_csf=3,
    )
    np.testing.assert_array_equal(speed[atlas > 0], np.float32(0.1))
    np.testing.assert_array_equal(speed[atlas <= 0], speed_before_merge[atlas <= 0])
    _assert_close(speed, expected_speed)

    result = ace_enhance_arrays(
        gray,
        white,
        multi_atlas_csf=atlas,
        intensity_threshold=126.9,
        tensor_alpha=1.6,
        gaussian_alpha=1.8,
        diffusion_iterations=25,
        use_updated_csf_in_ace=3,
    )
    _assert_close(result.enhanced_gray_matter, _array(REFERENCE / "enhanced_gm_java.nii.gz"))
    np.testing.assert_array_equal(result.skeleton, _array(REFERENCE / "csf_skeleton_java.nii.gz"))
    np.testing.assert_array_equal(result.thinned_skeleton, _array(REFERENCE / "thinned_skeleton_java.nii.gz"))
