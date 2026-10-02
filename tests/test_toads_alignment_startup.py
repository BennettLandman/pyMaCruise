from __future__ import annotations

import numpy as np
import nibabel as nib
import pytest

from pymacruise.toads_alignment_startup import (
    prepare_toads_alignment_startup,
    register_toads_alignment_startup,
)


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


def test_startup_memberships_feed_selected_exact_rigid_atlas_pyramid() -> None:
    t1, mask, topology, priors = _synthetic_inputs()
    startup = prepare_toads_alignment_startup(
        t1, mask, topology, priors, (1.0, 1.0, 1.0)
    )
    xx, yy, zz = np.indices((15, 13, 11), dtype=np.float32)
    centers = ((7, 6, 5), (5, 5, 4), (9, 7, 6), (6, 8, 5))
    atlas_shapes = np.empty((15, 13, 11, 4), dtype=np.float32)
    for klass, (cx, cy, cz) in enumerate(centers):
        distance = (xx - cx) ** 2 + 1.3 * (yy - cy) ** 2 + 0.8 * (zz - cz) ** 2
        atlas_shapes[..., klass] = np.exp(-distance / np.float32(18.0))
    atlas_shapes[..., 0] = np.float32(1.0)
    memberships_before = startup.initial_memberships.copy()

    result = register_toads_alignment_startup(
        startup,
        atlas_shapes,
        levels=2,
        max_iterations=1,
    )

    assert result.startup is startup
    assert len(result.registration.levels) == 1
    assert result.registration.levels[0].scale == 2
    assert result.registration.transform.shape == (6,)
    assert np.isfinite(result.registration.transform).all()
    setup = result.setup
    assert setup.memberships is startup.initial_memberships
    assert setup.best_candidates is None
    np.testing.assert_array_equal(setup.initial_zero_transform, np.zeros(6, dtype=np.float32))
    assert setup.transform is startup.transform
    assert setup.transform_type == "rigid"
    assert setup.approximate_class_count == priors.shape[1]
    assert setup.use_approximate_memberships is False
    assert setup.max_iterations == 1
    assert setup.levels == 2
    assert setup.subsample_factor == 2
    assert setup.signed_image_resolutions_mm == startup.geometry.signed_image_resolutions_mm
    assert setup.signed_image_resolutions_mm == (-1.0, -1.0, 1.0)
    assert setup.template_resolutions_mm == startup.template_resolution_mm
    assert setup.image_dimensions_xyz == startup.initial_memberships.shape[:3]
    assert setup.template_dimensions_xyz == (15, 13, 11)
    assert setup.mipav_image_orientation == 0
    assert setup.mipav_axis_orientation_codes == (2, 3, 5)
    assert setup.orientation_signs == (-1.0, -1.0, 1.0)
    assert setup.image_center_xyz == tuple(size / 2.0 for size in startup.initial_memberships.shape[:3])
    assert setup.template_center_xyz == (7.5, 6.5, 5.5)
    assert setup.stage_order == (
        "setImageInfo",
        "initShapeRegistration",
        "setTransform",
        "updateShapeRegistration",
        "registerShapesPyramid",
    )
    np.testing.assert_array_equal(startup.initial_memberships, memberships_before)


def test_registration_rejects_shape_class_order_mismatch() -> None:
    t1, mask, topology, priors = _synthetic_inputs()
    startup = prepare_toads_alignment_startup(
        t1, mask, topology, priors, (1.0, 1.0, 1.0)
    )
    with pytest.raises(ValueError, match="startup class order"):
        register_toads_alignment_startup(
            startup, np.ones((9, 8, 7, 3), dtype=np.float32), levels=2
        )
