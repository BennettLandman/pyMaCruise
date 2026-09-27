from __future__ import annotations

from pathlib import Path

import numpy as np
import nibabel as nib
import pytest

from pymacruise.tgdm import (
    TGDMProfileSettings,
    TGDMBaselineParameters,
    tgdm_evolve,
    tgdm_force,
    tgdm_geometry,
    tgdm_initial_membership_level_set,
    tgdm_inner_pressure_field,
    tgdm_central_pressure_field,
    tgdm_outer_pressure_field,
    tgdm_simple_point_618,
    tgdm_tissue_sign_field,
    tgdm_nested_arrays,
    tgdm_nested_images,
    _ras_vector_components,
)


def test_tgdm_geometry_matches_quadratic_derivatives_and_curvature_numerators() -> None:
    i, j, k = np.indices((8, 8, 8), dtype=np.float32)
    field = i * i + 2 * j * j + 3 * k * k

    geometry = tgdm_geometry(field)

    np.testing.assert_array_equal(geometry.gradient[3, 3, 3], [6.0, 12.0, 18.0])
    assert geometry.mean_curvature_numerator[3, 3, 3] == 3456.0
    assert geometry.gaussian_curvature_numerator[3, 3, 3] == 5184.0
    assert geometry.backward[3, 3, 3, 0] == 5.0
    assert geometry.forward[3, 3, 3, 0] == 7.0


def test_tgdm_flat_front_has_zero_curvature_and_expected_pressure_upwind() -> None:
    ramp = np.broadcast_to(np.arange(8, dtype=np.float32)[:, None, None], (8, 7, 6))
    geometry = tgdm_geometry(ramp)

    assert np.all(geometry.mean_curvature_numerator == 0)
    assert np.all(geometry.gaussian_curvature_numerator == 0)
    force = tgdm_force(
        geometry,
        profile="central",
        curvature_force=-0.15,
        pressure_force=1.5,
        region_force=np.ones(ramp.shape, dtype=np.float32),
    )
    np.testing.assert_allclose(force[3:7, 3, 3], 1.5)
    np.testing.assert_allclose(force[0, 3, 3], 0.0)


def test_tgdm_advection_uses_source_upwind_direction_and_external_weight() -> None:
    ramp = np.broadcast_to(np.arange(6, dtype=np.float32)[:, None, None], (6, 5, 4))
    geometry = tgdm_geometry(ramp)
    forward_advection = np.zeros((*ramp.shape, 3), dtype=np.float32)
    forward_advection[..., 0] = 1
    backward_advection = -forward_advection

    positive = tgdm_force(
        geometry,
        profile="central",
        curvature_force=0,
        pressure_force=0,
        external_force=2,
        advection=forward_advection,
    )
    negative = tgdm_force(
        geometry,
        profile="central",
        curvature_force=0,
        pressure_force=0,
        external_force=2,
        advection=backward_advection,
    )
    assert positive[3, 2, 2] == 2.0
    assert negative[3, 2, 2] == -2.0


def test_tgdm_pressure_upwind_selection_uses_signed_weighted_region_force() -> None:
    field = np.zeros((5, 5, 5), dtype=np.float32)
    field[1, 2, 2] = 1
    field[2, 2, 2] = 2
    field[3, 2, 2] = 4
    geometry = tgdm_geometry(field)

    positive = tgdm_force(
        geometry,
        profile="central",
        curvature_force=0,
        pressure_force=1.5,
        region_force=np.ones_like(field),
    )
    negative = tgdm_force(
        geometry,
        profile="central",
        curvature_force=0,
        pressure_force=-1.5,
        region_force=np.ones_like(field),
    )

    assert positive[2, 2, 2] == 1.5
    assert negative[2, 2, 2] == -3.0


def test_tgdm_rejects_invalid_geometry_and_force_fields() -> None:
    with pytest.raises(ValueError, match="3D"):
        tgdm_geometry(np.ones((5, 5)))
    with pytest.raises(ValueError, match="finite"):
        tgdm_geometry(np.full((4, 4, 4), np.inf))
    geometry = tgdm_geometry(np.zeros((4, 4, 4), dtype=np.float32))
    with pytest.raises(ValueError, match="region_force"):
        tgdm_force(
            geometry,
            profile="inner",
            curvature_force=1,
            pressure_force=1,
            region_force=np.ones((3, 3, 3)),
        )


def test_tgdm_geodesic_simple_point_accepts_one_interior_face_component() -> None:
    neighborhood = np.zeros((3, 3, 3), dtype=bool)
    neighborhood[0, 1, 1] = True

    assert tgdm_simple_point_618(neighborhood)


def test_tgdm_geodesic_simple_point_rejects_separate_interior_face_components() -> None:
    neighborhood = np.zeros((3, 3, 3), dtype=bool)
    neighborhood[0, 1, 1] = True
    neighborhood[2, 1, 1] = True

    assert not tgdm_simple_point_618(neighborhood)


def test_tgdm_geodesic_simple_point_rejects_missing_or_multiple_components() -> None:
    assert not tgdm_simple_point_618(np.zeros((3, 3, 3), dtype=bool))
    assert not tgdm_simple_point_618(np.ones((3, 3, 3), dtype=bool))


def test_tgdm_evolution_reinitializes_without_changing_a_planar_zero_force_front() -> None:
    i, _, _ = np.indices((9, 8, 7), dtype=np.float32)
    level_set = i - 4
    settings = TGDMProfileSettings(4.8, 3.6, 6.0, 20.0, 1, 2, 0.1, 0.045)

    result = tgdm_evolve(
        level_set,
        profile="inner",
        curvature_force=0,
        pressure_force=0,
        outer_iterations=1,
        settings=settings,
    )

    np.testing.assert_allclose(result.level_set, level_set, atol=1e-6)
    assert result.completed_steps == 2
    assert result.completed_outer_passes == 1
    assert result.landmine_restarts == 0


def test_tgdm_evolution_applies_inner_and_outer_surface_bounds() -> None:
    i, _, _ = np.indices((9, 7, 7), dtype=np.float32)
    level_set = i - 4
    inner_bound = level_set - 1
    outer_bound = level_set + 1
    settings = TGDMProfileSettings(4.8, 3.6, 6.0, 20.0, 0, 0, 0.1, 0.045)

    result = tgdm_evolve(
        level_set,
        profile="inner",
        curvature_force=0,
        pressure_force=0,
        inner_bound=inner_bound,
        outer_bound=outer_bound,
        outer_iterations=1,
        settings=settings,
    )

    # The inner bound is tighter on the negative side; the outer bound is
    # tighter on the positive side. Zero force leaves this nested interval.
    assert np.all(result.level_set <= outer_bound + 1e-6)
    assert np.all(result.level_set >= inner_bound - 1e-6)


def test_tgdm_evolution_moves_planar_front_with_signed_pressure() -> None:
    i, _, _ = np.indices((9, 8, 7), dtype=np.float32)
    level_set = i - 4
    settings = TGDMProfileSettings(4.8, 3.6, 6.0, 20.0, 0, 1, 0.1, 0.1)

    result = tgdm_evolve(
        level_set,
        profile="inner",
        curvature_force=0,
        pressure_force=1,
        outer_iterations=1,
        settings=settings,
    )

    assert result.level_set[4, 3, 3] == pytest.approx(-0.1, abs=1e-5)
    assert result.level_set[5, 3, 3] == pytest.approx(0.9, abs=1e-5)
    assert result.completed_steps == 1


def test_tgdm_initial_membership_level_set_uses_iso_minus_membership_sign() -> None:
    i, _, _ = np.indices((11, 7, 7), dtype=np.float32)
    membership = i * 40

    level_set = tgdm_initial_membership_level_set(membership, iso_level=160)

    assert level_set[3, 3, 3] > 0
    assert level_set[4, 3, 3] == pytest.approx(0, abs=1e-6)
    assert level_set[7, 3, 3] < 0
    assert np.max(np.abs(level_set)) <= 8


def test_tgdm_tissue_sign_field_matches_source_quantization_and_branches() -> None:
    wm = np.array([[[255, 0, 0, 128]]], dtype=np.float32)
    gm = np.array([[[0, 255, 0, 128]]], dtype=np.float32)

    result = tgdm_tissue_sign_field(wm, gm, white_matter_iso=0.1, gray_matter_iso=0.4)

    np.testing.assert_array_equal(result, np.array([[[254, 127, 0, 183]]], dtype=np.uint8))


def test_tgdm_surface_pressure_maps_match_source_scaling_and_skeleton_rule() -> None:
    wm = np.array([[[0, 25.5, 255]]], dtype=np.float32)
    np.testing.assert_allclose(tgdm_inner_pressure_field(wm), [[[-1, 0, 9]]])
    np.testing.assert_allclose(tgdm_central_pressure_field(np.array([[[0, 127, 254]]])), [[[-1, 0, 1]]])

    central = np.array([[[-1, 2, 5, 2]]], dtype=np.float32)
    gray = np.array([[[0, 101.6, 255, 0]]], dtype=np.float32)
    skeleton = np.array([[[0, 0, 0, 1]]], dtype=np.uint8)
    field = tgdm_outer_pressure_field(central, gray, skeleton, gray_matter_iso=0.4)
    np.testing.assert_allclose(field, [[[-1, 0, 35 / 34, -1.5]]])


def test_tgdm_nested_array_adapter_routes_tissues_and_supports_zero_iteration_smoke_case() -> None:
    x, y, z = np.indices((11, 11, 11))
    radius = np.sqrt((x - 5) ** 2 + (y - 5) ** 2 + (z - 5) ** 2)
    membership = np.where(radius <= 3, 255, 0).astype(np.float32)
    wm = np.where(radius <= 2, 255, 0).astype(np.float32)
    gm = np.where((radius > 2) & (radius <= 4), 150, 0).astype(np.float32)
    ace_gm = np.where((radius > 2) & (radius <= 4), 220, 0).astype(np.float32)
    parameters = TGDMBaselineParameters(
        initial_iterations=0,
        inner_iterations=0,
        central_shrink_iterations=0,
        central_iterations=0,
        outer_iterations=0,
    )

    result = tgdm_nested_arrays(membership, wm, gm, ace_gm, parameters=parameters)

    assert result.initial_level_set.shape == membership.shape
    assert result.initial_level_set[5, 5, 5] < 0
    assert result.tissue_sign.dtype == np.uint8
    assert result.tissue_sign[5, 5, 5] == 254
    for stage in (
        result.smoothed_initial,
        result.initial_inner,
        result.inner,
        result.central,
        result.central_shrink,
        result.outer,
    ):
        np.testing.assert_array_equal(stage.level_set, result.initial_level_set)
        assert stage.completed_steps == 0


def test_tgdm_image_adapter_reorients_warns_and_returns_ras_outputs() -> None:
    shape = (9, 8, 7)
    i, j, k = np.indices(shape)
    membership = np.where((i - 4) ** 2 + (j - 4) ** 2 + (k - 3) ** 2 <= 8, 255, 0).astype(np.float32)
    wm = np.where(membership > 0, 255, 0).astype(np.float32)
    gm = np.where(membership > 0, 100, 0).astype(np.float32)
    ace_gm = np.where(membership > 0, 220, 0).astype(np.float32)
    ras = np.eye(4)
    flipped = np.diag([-1.0, 1.0, 1.0, 1.0])
    flipped[0, 3] = shape[0] - 1
    params = TGDMBaselineParameters(
        initial_iterations=0,
        inner_iterations=0,
        central_shrink_iterations=0,
        central_iterations=0,
        outer_iterations=0,
    )

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        result = tgdm_nested_images(
            nib.Nifti1Image(membership, ras),
            nib.Nifti1Image(wm[::-1].copy(), flipped),
            nib.Nifti1Image(gm, ras),
            nib.Nifti1Image(ace_gm, ras),
            parameters=params,
        )

    for image in (
        result.initial_level_set,
        result.tissue_sign,
        result.initial_inner,
        result.inner,
        result.central,
        result.outer,
    ):
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert np.allclose(image.affine, ras)
    np.testing.assert_array_equal(result.execution.tissue_sign[4, 4, 3], 254)


def test_tgdm_image_adapter_rejects_mismatched_ras_grids() -> None:
    data = np.zeros((5, 5, 5), dtype=np.float32)
    reference = nib.Nifti1Image(data, np.eye(4))
    shifted = np.eye(4)
    shifted[0, 3] = 1
    params = TGDMBaselineParameters(
        initial_iterations=0,
        inner_iterations=0,
        central_shrink_iterations=0,
        central_iterations=0,
        outer_iterations=0,
    )

    with pytest.raises(ValueError, match="match the initial-membership RAS grid"):
        tgdm_nested_images(
            reference,
            nib.Nifti1Image(data, shifted),
            reference,
            reference,
            parameters=params,
        )


def test_tgdm_vector_components_follow_non_ras_axis_reorientation() -> None:
    shape = (5, 4, 3)
    affine = np.diag([-1.0, 1.0, 1.0, 1.0])
    affine[0, 3] = shape[0] - 1
    vector = np.zeros((*shape, 3), dtype=np.float32)
    vector[..., 0] = 1

    image = nib.Nifti1Image(vector, affine)
    ras = nib.as_closest_canonical(image)
    converted = _ras_vector_components(affine, np.asarray(ras.dataobj))

    assert nib.aff2axcodes(ras.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(converted[..., 0], -1)
    np.testing.assert_array_equal(converted[..., 1:], 0)


def test_tgdm_nested_synthetic_phantom_runs_every_evolution_stage() -> None:
    data_root = Path(__file__).parent / "data" / "synthetic" / "tgdm"
    images = {
        path.stem.removesuffix(".nii"): nib.load(str(path))
        for path in data_root.glob("*.nii.gz")
    }
    parameters = TGDMBaselineParameters(
        initial_iterations=1,
        inner_iterations=1,
        central_shrink_iterations=1,
        central_iterations=1,
        outer_iterations=1,
    )

    result = tgdm_nested_images(
        images["initial_membership"],
        images["white_matter"],
        images["gray_matter"],
        images["ace_gray_matter"],
        ace_skeleton=images["ace_skeleton"],
        parameters=parameters,
    )

    assert result.outer.shape == (25, 25, 25)
    for image in (result.initial_inner, result.inner, result.central, result.outer):
        assert np.isfinite(np.asarray(image.dataobj)).all()
    assert result.execution.initial_inner.completed_steps == 10
    assert result.execution.central.completed_steps == 10
    assert result.execution.central_shrink.completed_steps == 10
    assert result.execution.outer.completed_steps == 30
