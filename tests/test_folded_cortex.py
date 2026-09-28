from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
from scipy import ndimage

from pymacruise.ace import ace_enhance_arrays
from pymacruise.tgdm import TGDMProfileSettings, tgdm_evolve, tgdm_initial_membership_level_set


FOLDED = Path(__file__).parent / "data" / "synthetic" / "folded_cortex"


def _load(name: str) -> nib.Nifti1Image:
    return nib.load(FOLDED / f"{name}.nii.gz")


def test_ace_runs_on_larger_folded_synthetic_cortex() -> None:
    gray = np.asarray(_load("gray_matter").dataobj, dtype=np.float32)
    white = np.asarray(_load("white_matter").dataobj, dtype=np.float32)
    atlas_csf = np.asarray(_load("multi_atlas_csf").dataobj)

    result = ace_enhance_arrays(gray, white, multi_atlas_csf=atlas_csf)

    assert result.enhanced_gray_matter.shape == (41, 41, 41)
    for output in (
        result.enhanced_gray_matter,
        result.skeleton,
        result.thinned_skeleton,
    ):
        assert np.isfinite(output).all()
    assert set(np.unique(result.skeleton)) <= {0, 254}
    assert set(np.unique(result.thinned_skeleton)) <= {0, 254}
    assert np.count_nonzero(result.skeleton) > 0
    assert np.count_nonzero(result.thinned_skeleton) > 0


def test_tgdm_narrow_band_preserves_connected_folded_cortex_topology() -> None:
    membership = np.asarray(_load("initial_membership").dataobj, dtype=np.float32)
    seed = tgdm_initial_membership_level_set(membership, iso_level=126.9, max_distance=8.0)
    settings = TGDMProfileSettings(
        max_distance=4.7,
        how_close=3.6,
        max_signed_distance=6.2,
        final_max_signed_distance=20.0,
        internal_iterations=1,
        final_internal_iterations=2,
        delta_t=0.1,
        final_delta_t=0.045,
    )

    result = tgdm_evolve(
        seed,
        profile="outer",
        curvature_force=-1.0,
        pressure_force=0.0,
        outer_iterations=1,
        settings=settings,
    )
    mask = result.level_set <= 0

    assert result.completed_steps == 2
    assert np.isfinite(result.level_set).all()
    assert mask.any()
    assert ndimage.label(mask, structure=ndimage.generate_binary_structure(3, 1))[1] == 1
    assert ndimage.label(~mask, structure=ndimage.generate_binary_structure(3, 3))[1] == 1
