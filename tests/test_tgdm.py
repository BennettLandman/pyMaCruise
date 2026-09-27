from __future__ import annotations

import numpy as np
import pytest

from pymacruise.tgdm import tgdm_force, tgdm_geometry


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
