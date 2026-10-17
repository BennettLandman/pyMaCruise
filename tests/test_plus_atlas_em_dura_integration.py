from __future__ import annotations

import nibabel as nib
import numpy as np

from pymacruise.plus_atlas_em_dura import (
    run_atlas_em_prelesion_dura,
    run_dura_from_atlas_em_products,
)


def test_synthetic_atlasem_products_feed_selected_prelesion_dura_stage() -> None:
    shape = (17, 17, 17)
    x, y, z = np.indices(shape, dtype=np.float32)
    radius = np.sqrt((x - 8) ** 2 + (y - 8) ** 2 + (z - 8) ** 2)
    mask = (radius < 7.0).astype(np.uint8)
    membership = np.zeros((*shape, 3), dtype=np.float32)
    membership[..., 0] = np.where(mask, 0.05, 0.0)
    membership[..., 1] = np.where((radius >= 3.0) & (radius < 7.0), 0.75, 0.01)
    membership[..., 2] = np.where(radius < 3.0, 0.85, 0.02)
    membership[mask == 0] = 0.0
    t1_values = np.where(mask, 100.0 - 4.0 * radius, 0.0).astype(np.float32)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    t1 = nib.Nifti1Image(t1_values, affine)
    membership_image = nib.Nifti1Image(membership, affine)
    mask_image = nib.Nifti1Image(mask, affine)

    result = run_atlas_em_prelesion_dura(
        t1,
        membership_image,
        mask_image,
        ("Background", "CorticalGM", "WhiteMatter"),
        max_iterations=1,
    )

    assert result.products.shape == (*shape, 3)
    assert result.products.get_filename() == "atlasem_dura_inputs.nii.gz"
    assert nib.aff2axcodes(result.products.affine) == ("R", "A", "S")
    assert result.segmentation.shape == shape
    assert np.any(result.products.dataobj[..., 0] > 0)
    assert np.any(result.products.dataobj[..., 1] > 0)
    assert result.dura.evolution.completed_outer_passes == 1
    assert nib.aff2axcodes(result.dura.stripped_t1.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(result.dura.brain_mask.affine) == ("R", "A", "S")
    assert np.isfinite(np.asarray(result.dura.stripped_t1.dataobj)).all()
    file_handoff = run_dura_from_atlas_em_products(
        t1, result.products, max_iterations=1
    )
    assert file_handoff.products.shape == result.products.shape
    assert file_handoff.dura.evolution.completed_outer_passes == 1
    np.testing.assert_array_equal(
        file_handoff.dura.brain_mask.dataobj,
        result.dura.brain_mask.dataobj,
    )
