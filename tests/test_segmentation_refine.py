from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.segmentation_refine import refine_segmentation_by_surfaces


REFINEMENT_FIXTURES = Path(__file__).parent / "data" / "synthetic" / "refinement"


def _load_inputs():
    return tuple(
        nib.load(REFINEMENT_FIXTURES / f"{name}.nii.gz")
        for name in ("segmentation", "inner_levelset", "outer_levelset")
    )


def test_source_thresholds_and_all_refinement_rules_on_committed_phantom() -> None:
    segmentation, inner, outer = _load_inputs()

    result = refine_segmentation_by_surfaces(segmentation, inner, outer)
    output = np.asarray(result.image.dataobj)

    assert result.threshold_voxels == {
        "inner_wm": 0.2,
        "outer_gm": 0.2,
        "central_wm": 0.2,
        "central_gm": 0.2,
        "outer_central": 0.2,
    }
    assert output[3, 3, 3] == 44
    assert output[2, 8, 2] == 44  # Tie resolves to the lower label.
    assert output[8, 3, 3] == 0
    assert output[5, 5, 5] == 120
    assert output[7, 7, 7] == 120
    assert output[0, 0, 10] == 102  # No local WM vote; remains unresolved.
    assert output[10, 10, 10] == 7  # Unrelated labels remain untouched.
    assert result.changed_voxels == {
        "inner_wm": 2,
        "outer_clear": 1,
        "central_fill": 2,
    }
    assert result.unresolved_voxels == {
        "inner_wm": 1,
        "outer_clear": 0,
        "central_fill": 0,
    }
    assert result.sweeps["inner_wm"] >= 1
    assert nib.aff2axcodes(result.image.affine) == ("R", "A", "S")
    assert result.image.get_data_dtype() == segmentation.get_data_dtype()


def test_central_rule_can_be_disabled_and_outputs_explain_candidate_masks() -> None:
    segmentation, inner, outer = _load_inputs()

    result = refine_segmentation_by_surfaces(
        segmentation, inner, outer, deal_with_central=False
    )
    output = np.asarray(result.image.dataobj)
    assert output[5, 5, 5] == 45
    assert output[7, 7, 7] == 0
    assert result.changed_voxels["central_fill"] == 0
    assert not result.changed_masks["central_fill"].any()
    assert result.candidate_masks["inner_wm"][3, 3, 3]
    assert result.candidate_masks["outer_clear"][8, 3, 3]
    assert result.candidate_masks["central_wm"].sum() == 0


def test_non_ras_refinement_inputs_warn_and_produce_ras_result() -> None:
    ras_inputs = _load_inputs()
    expected = refine_segmentation_by_surfaces(*ras_inputs)
    non_ras = []
    for image in ras_inputs:
        affine = image.affine.copy()
        data = np.asarray(image.dataobj)[::-1]
        affine[:3, 0] *= -1
        affine[0, 3] += (image.shape[0] - 1) * image.affine[0, 0]
        non_ras.append(nib.Nifti1Image(data, affine))

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        actual = refine_segmentation_by_surfaces(*non_ras)
    np.testing.assert_array_equal(actual.image.dataobj, expected.image.dataobj)
    assert nib.aff2axcodes(actual.image.affine) == ("R", "A", "S")


def test_refinement_rejects_misaligned_and_noninteger_label_images() -> None:
    segmentation, inner, outer = _load_inputs()
    shifted_affine = outer.affine.copy()
    shifted_affine[0, 3] += 1
    shifted = nib.Nifti1Image(np.asarray(outer.dataobj), shifted_affine)
    with pytest.raises(ValueError, match="same RAS grid"):
        refine_segmentation_by_surfaces(segmentation, inner, shifted)

    fractional = nib.Nifti1Image(
        np.asarray(segmentation.dataobj, dtype=np.float32) + 0.25,
        segmentation.affine,
    )
    with pytest.raises(ValueError, match="integer labels"):
        refine_segmentation_by_surfaces(fractional, inner, outer)


def test_refinement_rejects_invalid_thresholds() -> None:
    segmentation, inner, outer = _load_inputs()
    with pytest.raises(ValueError, match="finite and nonnegative"):
        refine_segmentation_by_surfaces(
            segmentation, inner, outer, inner_coefficient_mm=-0.1
        )
    with pytest.raises(ValueError, match="radius must be nonnegative"):
        refine_segmentation_by_surfaces(
            segmentation, inner, outer, neighborhood_radius=-1
        )


def test_surface_coefficient_uses_largest_voxel_spacing() -> None:
    shape = (7, 7, 7)
    affine = np.diag([0.5, 0.5, 0.5, 1.0])
    labels = np.zeros(shape, dtype=np.uint16)
    labels[3, 3, 3] = 101
    labels[2:5, 2:5, 2:5] = 44
    labels[3, 3, 3] = 101
    inner = np.ones(shape, dtype=np.float32)
    inner[3, 3, 3] = -0.3
    outer = np.zeros(shape, dtype=np.float32)
    images = tuple(nib.Nifti1Image(data, affine) for data in (labels, inner, outer))

    below = refine_segmentation_by_surfaces(*images)
    assert below.threshold_voxels["inner_wm"] == pytest.approx(0.4)
    assert np.asarray(below.image.dataobj)[3, 3, 3] == 101

    inner[3, 3, 3] = -0.5
    above = refine_segmentation_by_surfaces(
        images[0], nib.Nifti1Image(inner, affine), images[2]
    )
    assert np.asarray(above.image.dataobj)[3, 3, 3] == 44
