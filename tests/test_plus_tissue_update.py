from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.boundary import combine_multi_atlas_toads_images
from pymacruise.plus_pipeline import run_plus_tissue_update_stage


FIXTURE = Path(__file__).parent / "data/synthetic/plus_tissue_update"
INPUT_NAMES = (
    "multi_atlas_gm",
    "multi_atlas_wm",
    "multi_atlas_brain_mask",
    "toads_gm",
    "toads_wm",
    "toads_mask",
)


def _inputs() -> list[nib.spatialimages.SpatialImage]:
    return [nib.load(FIXTURE / f"{name}.nii.gz") for name in INPUT_NAMES]


def test_plus_tissue_update_routes_six_inputs_to_selected_combiner() -> None:
    inputs = _inputs()
    actual = run_plus_tissue_update_stage(*inputs)
    expected = combine_multi_atlas_toads_images(*inputs, gm_confident_level=0.001)

    for output, direct in zip(
        (actual.gm, actual.wm, actual.mask), expected, strict=True
    ):
        np.testing.assert_array_equal(output.dataobj, direct.dataobj)
        np.testing.assert_array_equal(output.affine, inputs[3].affine)
        assert nib.aff2axcodes(output.affine) == ("R", "A", "S")
    assert [image.get_data_dtype() for image in (actual.gm, actual.wm, actual.mask)] == [
        np.dtype(np.float32),
        np.dtype(np.float32),
        np.dtype(np.int32),
    ]


def test_plus_tissue_update_preserves_selected_priority_threshold_and_mask_rules() -> None:
    actual = run_plus_tissue_update_stage(*_inputs())
    gm = np.asarray(actual.gm.dataobj)
    wm = np.asarray(actual.wm.dataobj)
    mask = np.asarray(actual.mask.dataobj)

    assert (gm[3, 3, 3], wm[3, 3, 3], mask[3, 3, 3]) == (0, 255, 255)
    assert (gm[1, 1, 1], wm[1, 1, 1], mask[1, 1, 1]) == (255, 0, 0)
    assert (gm[1, 2, 2], wm[1, 2, 2], mask[1, 2, 2]) == pytest.approx((0.255, 0.6, 7))
    assert (gm[0, 0, 0], wm[0, 0, 0], mask[0, 0, 0]) == (0, 0, 0)


def test_plus_tissue_update_rejects_physical_grid_mismatch() -> None:
    inputs = _inputs()
    shifted = inputs[-1].affine.copy()
    shifted[0, 3] += 1.0
    inputs[-1] = nib.Nifti1Image(np.asanyarray(inputs[-1].dataobj), shifted)
    with pytest.raises(ValueError, match="same physical grid"):
        run_plus_tissue_update_stage(*inputs)
