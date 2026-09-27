from __future__ import annotations

from pathlib import Path

import numpy as np
import nibabel as nib
import pytest
from scipy import ndimage

from pymacruise.ace import (
    _central_gradient_mirror,
    _distance_and_characteristic_fmm,
    _principal_tensor_vector,
    ace_anisotropic_diffusion,
    ace_enhance_arrays,
    ace_enhance,
    ace_close_wm,
    ace_gaussian_blur,
    ace_seed_fields,
    _signed_distance_fmm,
    ace_update_speed,
    clean_ace_components,
    narrow_band_ace,
)
from pymacruise.skeleton import thin_ace, thin_ace_image

ACE_PHANTOM = Path(__file__).parent / "data" / "synthetic" / "ace"


def test_deriche_gaussian_is_finite_repeatable_and_linear() -> None:
    shape = (7, 8, 9)
    first = np.zeros(shape, dtype=np.float32)
    first[3, 3, 4] = 1
    second = np.fromfunction(lambda i, j, k: i + 2 * j - k, shape).astype(np.float32)

    blurred = ace_gaussian_blur(first, 1.8)
    repeated = ace_gaussian_blur(first, 1.8)
    combined = ace_gaussian_blur(first + second, 1.8)

    assert blurred.dtype == np.float32
    assert np.isfinite(blurred).all()
    np.testing.assert_array_equal(blurred, repeated)
    np.testing.assert_allclose(combined, blurred + ace_gaussian_blur(second, 1.8), atol=2e-6)
    with pytest.raises(ValueError, match="positive"):
        ace_gaussian_blur(first, 0)


def test_ace_closing_uses_repeated_eighteen_neighbor_binary_morphology() -> None:
    data = np.zeros((13, 13, 13), dtype=np.float32)
    data[3:10, 3:10, 3:10] = 255
    data[6, 6, 6] = 0

    closed = ace_close_wm(data, iterations=1)
    assert closed.dtype == np.uint8
    assert set(np.unique(closed)) <= {0, 1}
    assert closed[6, 6, 6] == 1
    assert closed[0, 0, 0] == 0
    np.testing.assert_array_equal(ace_close_wm(data, iterations=0), (data != 0).astype(np.uint8))


def test_ace_cleanup_removes_small_interior_components_but_keeps_border_only() -> None:
    data = np.zeros((12, 12, 12), dtype=np.uint8)
    data[2:4, 2:4, 2:4] = 254  # 8 voxels: below the threshold.
    data[5:10, 5:10, 5:10] = 254  # 125 voxels: retained.
    data[0, 0:3, 0:3] = 254  # Border-only component has no source seed voxel.

    cleaned = clean_ace_components(data, minimum_size_exclusive=60)

    assert not np.any(cleaned[2:4, 2:4, 2:4])
    assert np.all(cleaned[5:10, 5:10, 5:10] == 254)
    assert np.all(cleaned[0, 0:3, 0:3] == 254)


def test_ace_helpers_reject_non_3d_arrays() -> None:
    with pytest.raises(ValueError, match="3D"):
        ace_close_wm(np.ones((4, 4)))
    with pytest.raises(ValueError, match="3D"):
        clean_ace_components(np.ones((4, 4)))


def test_ace_seed_fields_initialize_hollow_white_matter_with_signed_distance() -> None:
    gray = np.zeros((15, 15, 15), dtype=np.float32)
    white = np.zeros_like(gray)
    white[4:11, 4:11, 4:11] = 220
    white[6:9, 6:9, 6:9] = 0
    phi = np.ones_like(gray)

    fields = ace_seed_fields(gray, white, phi)

    assert fields.initial_csf_speed[7, 7, 7] == 1.0
    assert fields.prohibited_region[7, 7, 7] == 0
    assert fields.prohibited_region[5, 5, 5] == 1
    assert fields.naive_csf[7, 7, 7] == 254
    assert fields.naive_csf[5, 5, 5] == 34


def test_ace_anisotropic_diffusion_is_stable_bounded_and_monotone() -> None:
    shape = (7, 8, 9)
    initial = np.zeros(shape, dtype=np.float32)
    initial[2:5, 3:6, 3:6] = 0.4
    guide = np.fromfunction(lambda i, j, k: 2 * i + j - k, shape).astype(np.float32)

    diffused = ace_anisotropic_diffusion(initial, guide, iterations=4)
    repeated = ace_anisotropic_diffusion(initial, guide, iterations=4)

    assert np.isfinite(diffused).all()
    assert diffused.dtype == np.float32
    assert diffused.min() >= 0.0
    assert diffused.max() <= 1.0
    assert np.all(diffused >= initial)
    np.testing.assert_array_equal(diffused, repeated)
    np.testing.assert_array_equal(
        ace_anisotropic_diffusion(initial, guide, iterations=0), initial
    )


def test_source_difference_edges_and_power_method_orientation() -> None:
    ramp = np.broadcast_to(np.arange(4, dtype=np.float32)[:, None, None], (4, 3, 2))
    gx, gy, gz = _central_gradient_mirror(ramp)
    np.testing.assert_array_equal(gx[:, 0, 0], [0.5, 1, 1, 0.5])
    assert np.all(gy == 0) and np.all(gz == 0)

    one = np.ones((2, 2, 2), dtype=np.float32)
    zero = np.zeros_like(one)
    vx, vy, vz = _principal_tensor_vector(9 * one, zero, zero, 4 * one, zero, one)
    np.testing.assert_allclose(vx, 1.0)
    np.testing.assert_allclose(vy, 0.0)
    np.testing.assert_allclose(vz, 0.0)


def test_updated_multi_atlas_csf_is_merged_after_diffusion() -> None:
    csf = np.zeros((5, 6, 7), dtype=np.float32)
    csf[2, 3, 4] = 0.25
    atlas = np.zeros_like(csf)
    atlas[1, 2, 3] = 1
    merged = csf.copy()
    merged[atlas > 0] = 1.0

    expected = ace_gaussian_blur(1.0 - 0.9 * merged, 1.8)
    actual = ace_update_speed(csf, multi_atlas_csf=atlas)
    np.testing.assert_array_equal(actual, expected)


def test_ace_signed_distance_matches_subvoxel_plane_and_caps_empty_interface() -> None:
    phi = np.broadcast_to(
        (np.arange(9, dtype=np.float32) - 3.25)[:, None, None], (9, 7, 6)
    ).copy()
    distance = _signed_distance_fmm(phi)
    np.testing.assert_allclose(distance, phi, atol=1e-6)
    assert distance[0, 0, 0] == -3.25
    assert distance[-1, 0, 0] == 4.75

    no_interface = _signed_distance_fmm(np.ones((3, 4, 5), dtype=np.float32), 12.0)
    np.testing.assert_array_equal(no_interface, 12.0)


def test_narrow_band_ace_updates_synchronously_and_reinitializes_to_four_voxels() -> None:
    phi = np.broadcast_to(
        (np.arange(10, dtype=np.float32) - 4.25)[:, None, None], (10, 8, 7)
    ).copy()

    evolved = narrow_band_ace(phi)

    assert np.isfinite(evolved).all()
    assert evolved.dtype == np.float32
    assert np.max(evolved) == 4.0
    assert np.min(evolved) >= -4.0
    # The planar front remains planar, with the source's positive-speed update
    # translating its subvoxel crossing by approximately 0.001 voxel.
    np.testing.assert_allclose(evolved[:, 3, 3], evolved[:, 4, 4], atol=1e-6)
    assert abs(evolved[4, 3, 3] - phi[4, 3, 3]) < 0.01


def test_ace_speed_weighted_fmm_scales_plane_distance_and_has_characteristics() -> None:
    phi = np.broadcast_to(
        (np.arange(8, dtype=np.float32) - 3.25)[:, None, None], (8, 6, 5)
    ).copy()
    speed = np.full(phi.shape, 0.5, dtype=np.float32)

    travel_time, characteristic = _distance_and_characteristic_fmm(phi, speed)

    np.testing.assert_allclose(travel_time[:, 2, 2], phi[:, 2, 2] * 2.0, atol=1e-6)
    assert np.isfinite(characteristic).all()
    assert np.all(characteristic >= 0)
    empty_time, empty_char = _distance_and_characteristic_fmm(
        np.ones((3, 4, 5), dtype=np.float32), speed=np.ones((3, 4, 5), dtype=np.float32)
    )
    assert not np.any(empty_time)
    assert not np.any(empty_char)


def test_ace_array_pipeline_is_deterministic_and_clears_volume_borders() -> None:
    gray = np.asarray(nib.load(ACE_PHANTOM / "gray_matter.nii.gz").dataobj)
    white = np.asarray(nib.load(ACE_PHANTOM / "white_matter.nii.gz").dataobj)
    atlas = np.asarray(nib.load(ACE_PHANTOM / "multi_atlas_csf.nii.gz").dataobj)

    result = ace_enhance_arrays(
        gray, white, multi_atlas_csf=atlas, diffusion_iterations=2
    )
    repeated = ace_enhance_arrays(
        gray, white, multi_atlas_csf=atlas, diffusion_iterations=2
    )

    assert result.enhanced_gray_matter.shape == gray.shape
    assert result.enhanced_gray_matter.dtype == np.float32
    assert result.skeleton.dtype == np.int32
    assert set(np.unique(result.skeleton)) <= {0, 254}
    assert np.isfinite(result.enhanced_gray_matter).all()
    assert np.all(result.enhanced_gray_matter[[0, -1], :, :] == 0)
    assert np.all(result.skeleton[:, [0, -1], :] == 0)
    assert np.count_nonzero(result.skeleton) > 0
    np.testing.assert_array_equal(result.enhanced_gray_matter, repeated.enhanced_gray_matter)
    np.testing.assert_array_equal(result.skeleton, repeated.skeleton)


def test_ace_nifti_api_warns_reorients_non_ras_and_returns_ras_outputs() -> None:
    gray = np.asarray(nib.load(ACE_PHANTOM / "gray_matter.nii.gz").dataobj)
    white = np.asarray(nib.load(ACE_PHANTOM / "white_matter.nii.gz").dataobj)
    shape = gray.shape
    affine = np.diag([-1.0, 1.0, 1.0, 1.0])
    affine[0, 3] = shape[0] - 1
    gray_image = nib.Nifti1Image(gray, affine)
    white_image = nib.Nifti1Image(white, affine)

    with pytest.warns(UserWarning, match="reorienting"):
        result = ace_enhance(gray_image, white_image, diffusion_iterations=1)

    assert nib.aff2axcodes(result.enhanced_gray_matter.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(result.skeleton.affine) == ("R", "A", "S")
    np.testing.assert_allclose(result.enhanced_gray_matter.affine, np.eye(4))
    np.testing.assert_allclose(result.skeleton.affine, np.eye(4))


def test_thin_ace_preserves_line_and_reduces_solid_while_clearing_border() -> None:
    line = np.zeros((15, 15, 15), dtype=np.uint8)
    line[2:13, 7, 7] = 255
    thin_line = thin_ace(line)
    np.testing.assert_array_equal(thin_line[2:13, 7, 7], 254)
    assert np.count_nonzero(thin_line) == np.count_nonzero(line)

    solid = np.zeros_like(line)
    solid[2:13, 2:13, 2:13] = 255
    thin_solid = thin_ace(solid)
    assert 0 < np.count_nonzero(thin_solid) < np.count_nonzero(solid)
    assert set(np.unique(thin_solid)) <= {0, 254}
    assert np.all(thin_solid[[0, -1], :, :] == 0)
    assert ndimage.label(
        thin_solid > 0, structure=ndimage.generate_binary_structure(3, 2)
    )[1] == 1


def test_thin_ace_nifti_wrapper_normalizes_to_ras() -> None:
    image = nib.load(ACE_PHANTOM / "multi_atlas_csf.nii.gz")
    flipped_affine = image.affine.copy()
    flipped_affine[0, 0] *= -1
    flipped_affine[0, 3] = image.shape[0] - 1
    non_ras = nib.Nifti1Image(np.asarray(image.dataobj), flipped_affine)

    with pytest.warns(UserWarning, match="reorienting"):
        result = thin_ace_image(non_ras)

    assert nib.aff2axcodes(result.skeleton.affine) == ("R", "A", "S")
