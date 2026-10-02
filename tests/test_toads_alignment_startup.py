from __future__ import annotations

import numpy as np
import nibabel as nib
import pytest

from pymacruise.toads_alignment_startup import prepare_toads_alignment_startup


def _synthetic_inputs():
    shape = (28, 27, 25)
    x, y, z = np.indices(shape, dtype=np.float32)
    mask_values = np.zeros(shape, dtype=np.uint8)
    mask_values[5:23, 4:23, 3:22] = 1
    values = np.zeros(shape, dtype=np.float32)
    inside = mask_values > 0
    values[inside] = (
        np.float32(20.0)
        + x[inside] * np.float32(0.7)
        + y[inside] * np.float32(0.31)
        + z[inside] * np.float32(0.19)
        + ((x[inside] + 2 * y[inside] + z[inside]) % 7) * np.float32(0.13)
    )
    affine = np.diag((1.0, 1.0, 1.0, 1.0))
    t1 = nib.Nifti1Image(values, affine)
    mask = nib.Nifti1Image(mask_values, affine)

    # Asymmetric, folded synthetic atlas support with source background label 1.
    atlas = np.ones((15, 13, 11), dtype=np.uint8)
    for i in range(2, 13):
        for j in range(2, 11):
            for k in range(2, 9):
                fold = (i + 2 * j + k) % 5
                if (i - 7) ** 2 + (j - 6) ** 2 + (k - 5) ** 2 < 24 + fold:
                    atlas[i, j, k] = 2 + ((i + j + k) % 3)
    priors = np.asarray([[0.10, 0.36, 0.63, 0.88]], dtype=np.float32)
    return t1, mask, atlas, priors


def test_startup_composes_selected_pre_registration_stages_in_order() -> None:
    t1, mask, topology, priors = _synthetic_inputs()

    result = prepare_toads_alignment_startup(
        t1, mask, topology, priors, (1.0, 1.0, 1.0)
    )

    assert result.stage_order == (
        "normalizeImages",
        "findCroppingBoundaries",
        "cropImages",
        "initialCentroids",
        "initialMemberships",
    )
    assert result.intensity.image.shape == t1.shape
    assert result.cropped_images.shape[0] == 1
    assert result.cropped_topology.shape == result.cropped_images.shape[1:]
    assert result.initial_memberships.shape == (*result.cropped_images.shape[1:], priors.shape[1])
    assert result.transform.shape == (6,)
    assert result.geometry.orientation_signs == (-1.0, -1.0, 1.0)
    assert result.crop_bounds[0] <= result.crop_bounds[1]
    assert result.crop_bounds[2] <= result.crop_bounds[3]
    assert result.crop_bounds[4] <= result.crop_bounds[5]

    # The exact initializer leaves crop padding and each terminal source plane zero.
    border = 10
    assert np.count_nonzero(result.initial_memberships[:border]) == 0
    assert np.count_nonzero(result.initial_memberships[:, :border]) == 0
    assert np.count_nonzero(result.initial_memberships[:, :, :border]) == 0
    assert np.count_nonzero(result.initial_memberships[-1]) == 0
    assert np.count_nonzero(result.initial_memberships[:, -1]) == 0
    assert np.count_nonzero(result.initial_memberships[:, :, -1]) == 0


def test_startup_reorients_non_ras_t1_and_mask_with_warning() -> None:
    t1, mask, topology, priors = _synthetic_inputs()
    flip = np.eye(4)
    flip[0, 0] = -1
    flip[1, 1] = -1
    flip[0, 3] = t1.shape[0] - 1
    flip[1, 3] = t1.shape[1] - 1
    non_ras_t1 = nib.Nifti1Image(np.asarray(t1.dataobj)[::-1, ::-1, :], t1.affine @ flip)
    non_ras_mask = nib.Nifti1Image(np.asarray(mask.dataobj)[::-1, ::-1, :], mask.affine @ flip)

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        result = prepare_toads_alignment_startup(
            non_ras_t1, non_ras_mask, topology, priors, (1.0, 1.0, 1.0)
        )

    assert nib.aff2axcodes(result.intensity.image.affine) == ("R", "A", "S")
    assert result.input_orientation_codes == ("L", "P", "S")


@pytest.mark.parametrize(
    "priors,template_resolution",
    [
        (np.ones((2, 3), dtype=np.float32), (1.0, 1.0, 1.0)),
        (np.ones((1, 3), dtype=np.float32), (1.0, 0.0, 1.0)),
    ],
)
def test_startup_rejects_invalid_atlas_contract(priors, template_resolution) -> None:
    t1, mask, topology, _ = _synthetic_inputs()
    with pytest.raises(ValueError):
        prepare_toads_alignment_startup(
            t1, mask, topology, priors, template_resolution
        )
