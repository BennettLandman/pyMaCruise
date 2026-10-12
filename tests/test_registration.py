from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.registration import (
    RegistrationConfig,
    apply_image_transform,
    apply_transforms_to_points_ras,
    apply_transforms_to_surface_ras,
    register_flair_to_t1,
    register_to_reference,
)
from pymacruise.surfaces import SurfaceMesh


REGISTRATION_FIXTURES = Path(__file__).parent / "data" / "synthetic" / "registration"


def _load(name: str):
    return nib.load(REGISTRATION_FIXTURES / f"{name}.nii.gz")


def _transform(name: str) -> Path:
    return REGISTRATION_FIXTURES / name


def test_known_forward_transform_resamples_intensity_and_labels_in_correct_direction() -> None:
    fixed = _load("fixed_image")
    moving = _load("moving_image")
    transform = _transform("known_moving_to_fixed.mat")

    warped_intensity = apply_image_transform(fixed, moving, [transform])
    warped_labels = apply_image_transform(
        _load("fixed_labels"),
        _load("moving_labels"),
        [transform],
        interpolation="label",
    )

    assert warped_intensity.shape == fixed.shape
    assert nib.aff2axcodes(warped_intensity.affine) == ("R", "A", "S")
    np.testing.assert_allclose(
        np.asarray(warped_intensity.dataobj)[2:-2, 2:-2, 2:-2],
        np.asarray(fixed.dataobj)[2:-2, 2:-2, 2:-2],
        atol=0.015,
        rtol=0,
    )
    np.testing.assert_array_equal(warped_labels.dataobj, _load("fixed_labels").dataobj)
    assert warped_labels.get_data_dtype() == np.uint8


def test_inverse_transform_restores_moving_grid_and_point_coordinates() -> None:
    moving = _load("moving_labels")
    fixed = _load("fixed_labels")
    inverse = _transform("known_fixed_to_moving.mat")
    warped = apply_image_transform(
        moving,
        fixed,
        [inverse],
        interpolation="label",
    )
    np.testing.assert_array_equal(warped.dataobj, moving.dataobj)

    points = np.array([[5.0, 8.0, 7.0], [7.0, 8.0, 7.0]])
    forward_points = apply_transforms_to_points_ras(
        points, [_transform("known_moving_to_fixed.mat")]
    )
    inverse_points = apply_transforms_to_points_ras(
        forward_points, [inverse]
    )
    np.testing.assert_allclose(forward_points, points + [2, 0, 0], atol=1e-8)
    np.testing.assert_allclose(inverse_points, points, atol=1e-8)


def test_surface_transform_maps_ras_vertices_and_preserves_mesh_topology() -> None:
    vertices = np.array(
        [[4.0, 7.0, 6.0], [6.0, 7.0, 6.0], [5.0, 9.0, 6.0], [5.0, 8.0, 8.0]]
    )
    mesh = SurfaceMesh(
        vertices,
        np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]]),
    )
    transformed = apply_transforms_to_surface_ras(
        mesh, [_transform("known_moving_to_fixed.mat")]
    )
    np.testing.assert_allclose(transformed.vertices, vertices + [2, 0, 0], atol=1e-8)
    np.testing.assert_array_equal(transformed.faces, mesh.faces)


def test_registration_mask_and_estimated_rigid_transform_recover_known_shift(tmp_path: Path) -> None:
    result = register_to_reference(
        _load("fixed_image"),
        _load("moving_image"),
        fixed_mask=_load("fixed_mask"),
        output_prefix=tmp_path / "registration" / "fixture_",
        config=RegistrationConfig(
            transform_type="Rigid",
            metric="meansquares",
            sampling_rate=1.0,
            iterations=(80, 40, 0),
            shrink_factors=(2, 1, 1),
            smoothing_sigmas=(1, 0, 0),
        ),
    )

    assert result.transform_type == "Rigid"
    assert result.forward_transforms and result.inverse_transforms
    assert all(Path(path).is_file() for path in result.forward_transforms)
    assert nib.aff2axcodes(result.warped_moving.affine) == ("R", "A", "S")
    peak = np.unravel_index(np.argmax(result.warped_moving.get_fdata()), result.warped_moving.shape)
    assert np.linalg.norm(np.asarray(peak) - [7, 8, 7]) <= 1.0
    point = apply_transforms_to_points_ras(
        np.array([[5.0, 8.0, 7.0]]), result.forward_transforms
    )
    np.testing.assert_allclose(point, [[7.0, 8.0, 7.0]], atol=1.0)


def test_flair_to_t1_adapter_is_rigid_and_recovers_synthetic_shift(tmp_path: Path) -> None:
    result = register_flair_to_t1(
        _load("fixed_image"),
        _load("moving_flair"),
        fixed_mask=_load("fixed_mask"),
        output_prefix=tmp_path / "flair_to_t1_",
        config=RegistrationConfig(
            transform_type="Rigid",
            metric="mattes",
            sampling_rate=1.0,
            iterations=(60, 40, 20, 0),
            shrink_factors=(6, 4, 2, 1),
            smoothing_sigmas=(3, 2, 1, 0),
        ),
    )

    point = apply_transforms_to_points_ras(
        np.array([[5.0, 8.0, 7.0]]), result.forward_transforms
    )
    np.testing.assert_allclose(point, [[7.0, 8.0, 7.0]], atol=1.0)
    assert result.transform_type == "Rigid"
    assert all(Path(path).is_file() for path in result.forward_transforms)
    assert all(Path(path).is_file() for path in result.inverse_transforms)
    assert nib.aff2axcodes(result.warped_moving.affine) == ("R", "A", "S")
    with pytest.raises(ValueError, match="must be rigid-only"):
        register_flair_to_t1(
            _load("fixed_image"),
            _load("moving_flair"),
            output_prefix=tmp_path / "affine_",
            config=RegistrationConfig(transform_type="Affine"),
        )


def test_rigid_registration_repeats_with_fixed_seed(tmp_path: Path) -> None:
    config = RegistrationConfig(
        transform_type="Rigid",
        metric="mattes",
        sampling_rate=0.25,
        iterations=(60, 40, 20, 0),
        shrink_factors=(6, 4, 2, 1),
        smoothing_sigmas=(3, 2, 1, 0),
        deterministic=True,
        random_seed=321,
    )
    first = register_flair_to_t1(
        _load("fixed_image"),
        _load("moving_flair"),
        output_prefix=tmp_path / "repeat_a_",
        config=config,
    )
    second = register_flair_to_t1(
        _load("fixed_image"),
        _load("moving_flair"),
        output_prefix=tmp_path / "repeat_b_",
        config=config,
    )

    np.testing.assert_array_equal(
        np.asarray(first.warped_moving.dataobj),
        np.asarray(second.warped_moving.dataobj),
    )
    assert Path(first.forward_transforms[0]).read_bytes() == Path(
        second.forward_transforms[0]
    ).read_bytes()


def test_non_ras_inputs_warn_and_registration_grid_mismatches_fail() -> None:
    fixed = _load("fixed_image")
    moving = _load("moving_image")
    affine = moving.affine.copy()
    data = np.asarray(moving.dataobj)[::-1]
    affine[0, 0] = -1
    affine[0, 3] = moving.shape[0] - 1
    non_ras = nib.Nifti1Image(data, affine)
    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        result = apply_image_transform(
            fixed, non_ras, [_transform("known_moving_to_fixed.mat")]
        )
    assert nib.aff2axcodes(result.affine) == ("R", "A", "S")

    bad_mask_affine = fixed.affine.copy()
    bad_mask_affine[0, 3] += 1
    bad_mask = nib.Nifti1Image(np.asarray(_load("fixed_mask").dataobj), bad_mask_affine)
    with pytest.raises(ValueError, match="share the fixed RAS grid"):
        register_to_reference(
            fixed,
            moving,
            fixed_mask=bad_mask,
            output_prefix=Path(".pytest-tmp") / "bad_mask_",
        )


def test_registration_rejects_missing_transforms_and_noninteger_labels() -> None:
    fixed = _load("fixed_image")
    moving = _load("moving_image")
    with pytest.raises(FileNotFoundError, match="transform file does not exist"):
        apply_image_transform(fixed, moving, [REGISTRATION_FIXTURES / "absent.mat"])

    fractional = nib.Nifti1Image(
        _load("moving_image").get_fdata() * 3.7, moving.affine
    )
    with pytest.raises(ValueError, match="integer-valued"):
        apply_image_transform(
            fixed,
            fractional,
            [_transform("known_moving_to_fixed.mat")],
            interpolation="label",
        )
