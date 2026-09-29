from __future__ import annotations

from dataclasses import replace

import nibabel as nib
import numpy as np
import pytest

from pymacruise.cropping import crop_image, find_crop_parameters, uncrop_image
from pymacruise.plus_pipeline import uncrop_plus_tgdm_outputs
from pymacruise.tgdm import TGDMBaselineParameters, tgdm_nested_images


def _result_and_crop():
    shape = (9, 9, 9)
    affine = np.diag([1.2, 1.4, 1.6, 1.0])
    i, j, k = np.indices(shape)
    roi = (i - 4) ** 2 + (j - 4) ** 2 + (k - 4) ** 2 <= 8
    membership = np.where(roi, 255, 0).astype(np.float32)
    parameters = find_crop_parameters(
        nib.Nifti1Image(roi.astype(np.uint8), affine), edge_width=1
    )

    def crop(data: np.ndarray) -> nib.Nifti1Image:
        return crop_image(nib.Nifti1Image(data, affine), parameters)

    initial = crop(membership)
    wm = crop(np.where(roi, 255, 0).astype(np.float32))
    gm = crop(np.where(roi, 100, 0).astype(np.float32))
    ace_gm = crop(np.where(roi, 220, 0).astype(np.float32))
    skeleton_data = np.zeros(shape, dtype=np.uint8)
    skeleton_data[4, 4, 4] = 255
    skeleton = crop(skeleton_data)
    gvf = nib.Nifti1Image(
        np.zeros((*parameters.cropped_shape, 3), dtype=np.float32), initial.affine
    )
    zero_steps = TGDMBaselineParameters(
        initial_iterations=0,
        inner_iterations=0,
        central_shrink_iterations=0,
        central_iterations=0,
        outer_iterations=0,
    )
    result = tgdm_nested_images(
        initial,
        wm,
        gm,
        ace_gm,
        ace_skeleton=skeleton,
        gvf=gvf,
        parameters=zero_steps,
    )
    return result, parameters


def test_plus_tgdm_outputs_uncrop_to_full_ras_grid_with_source_max_rule() -> None:
    result, parameters = _result_and_crop()
    actual = uncrop_plus_tgdm_outputs(result, parameters)

    for name in (
        "initial_level_set",
        "tissue_sign",
        "initial_inner",
        "inner",
        "central",
        "outer",
    ):
        image = getattr(actual, name)
        expected = uncrop_image(
            getattr(result, name), parameters, legacy_exclusive_max=True
        )
        np.testing.assert_array_equal(image.dataobj, expected.dataobj)
        np.testing.assert_array_equal(image.affine, parameters.original_affine)
        assert image.shape == parameters.original_shape
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")

    cropped = np.asarray(result.inner.dataobj)
    uncropped = np.asarray(actual.inner.dataobj)
    assert uncropped[6, 6, 6] == cropped[5, 5, 5]
    assert uncropped[7, 7, 7] == parameters.background


def test_plus_tgdm_uncrop_rejects_a_product_on_another_crop_grid() -> None:
    result, parameters = _result_and_crop()
    shifted = result.inner.affine.copy()
    shifted[0, 3] += 1.0
    bad_inner = nib.Nifti1Image(np.asarray(result.inner.dataobj), shifted)
    bad_result = replace(result, inner=bad_inner)
    with pytest.raises(ValueError, match="inner must match the cropped RAS grid"):
        uncrop_plus_tgdm_outputs(bad_result, parameters)
