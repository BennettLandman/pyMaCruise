from __future__ import annotations

import numpy as np
import pytest

from pymacruise.tgdm import (
    TGDMProfileSettings,
    tgdm_evolve,
    tgdm_force,
    tgdm_geometry,
    tgdm_simple_point_618,
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
