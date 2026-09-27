from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.surfaces import (
    SurfaceMesh,
    extract_levelset_surface,
    read_legacy_vtk,
    surface_vertex_normals,
    transform_surface_ras,
    write_legacy_vtk,
)
from pymacruise.thickness import closest_point_thickness, levelset_subtraction_thickness


THICKNESS_FIXTURES = (
    Path(__file__).parent / "data" / "synthetic" / "surfaces" / "thickness"
)


def _load_surface(name: str) -> SurfaceMesh:
    image = nib.load(THICKNESS_FIXTURES / f"{name}.nii.gz")
    return extract_levelset_surface(image)


def test_levelset_subtraction_matches_constant_known_sphere_shell() -> None:
    inner = nib.load(THICKNESS_FIXTURES / "inner.nii.gz")
    outer = nib.load(THICKNESS_FIXTURES / "outer.nii.gz")

    result = levelset_subtraction_thickness(inner, outer)

    values = np.asarray(result.image.dataobj)
    assert values.shape == inner.shape
    assert nib.aff2axcodes(result.image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(values, 2.0, atol=5e-7)
    assert result.crossing_voxel_count > 0
    assert result.mean_mm == 2.0
    assert result.standard_deviation_mm == 0.0
    assert result.minimum_mm == result.maximum_mm == 2.0


def test_levelset_subtraction_warns_on_non_ras_and_rejects_grid_mismatch() -> None:
    inner = nib.load(THICKNESS_FIXTURES / "inner.nii.gz")
    data = np.asarray(inner.dataobj)[::-1]
    affine = inner.affine.copy()
    affine[0, 0] *= -1
    affine[0, 3] = (inner.shape[0] - 1) * 0.5
    non_ras_inner = nib.Nifti1Image(data, affine)
    outer = nib.load(THICKNESS_FIXTURES / "outer.nii.gz")
    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        result = levelset_subtraction_thickness(non_ras_inner, outer)
    assert nib.aff2axcodes(result.image.affine) == ("R", "A", "S")

    shifted_affine = outer.affine.copy()
    shifted_affine[0, 3] += 1.0
    shifted = nib.Nifti1Image(np.asarray(outer.dataobj), shifted_affine)
    with pytest.raises(ValueError, match="same RAS grid"):
        levelset_subtraction_thickness(inner, shifted)


def test_closest_point_thickness_reports_bidirectional_triangle_distances() -> None:
    inner = _load_surface("inner")
    outer = _load_surface("outer")
    central = _load_surface("central")

    result = closest_point_thickness(inner, outer, central=central)

    inner_data = result.inner_surface.point_data
    assert set(inner_data) == {
        "distance_inner_to_outer_mm",
        "distance_outer_to_inner_mm",
        "mean_thickness_mm",
    }
    assert np.mean(inner_data["distance_inner_to_outer_mm"]) == pytest.approx(2.0, abs=0.08)
    assert np.mean(inner_data["distance_outer_to_inner_mm"]) == pytest.approx(2.0, abs=0.08)
    assert result.mean_mm == pytest.approx(2.0, abs=0.08)
    assert np.max(np.abs(inner_data["mean_thickness_mm"] - 2.0)) < 0.15
    assert result.standard_deviation_mm < 0.12
    assert result.nonzero_vertex_count == len(inner.vertices)

    central_data = result.central_surface.point_data
    assert set(central_data) == {
        "distance_inner_to_outer_mm",
        "distance_central_to_inner_mm",
        "distance_to_inner_vertex_mm",
        "nearest_inner_vertex_index",
    }
    assert np.mean(central_data["distance_inner_to_outer_mm"]) == pytest.approx(2.0, abs=0.12)
    assert np.mean(central_data["distance_central_to_inner_mm"]) == pytest.approx(1.0, abs=0.08)
    assert np.all(central_data["nearest_inner_vertex_index"] >= 0)
    assert np.all(central_data["nearest_inner_vertex_index"] < len(inner.vertices))


def test_surface_normals_are_unit_outward_ras_vectors() -> None:
    mesh = _load_surface("inner")
    normals = surface_vertex_normals(mesh)
    center = np.mean(mesh.vertices, axis=0)

    np.testing.assert_allclose(np.linalg.norm(normals, axis=1), 1.0, atol=1e-12)
    assert np.mean(np.einsum("ij,ij->i", normals, mesh.vertices - center)) > 0


def test_full_space_surface_transform_uses_ras_mm_and_preserves_orientation() -> None:
    mesh = _load_surface("inner")
    transform = np.eye(4)
    transform[:3, 3] = [12.0, -4.0, 2.5]
    moved = transform_surface_ras(mesh, transform)
    np.testing.assert_allclose(moved.vertices, mesh.vertices + [12.0, -4.0, 2.5])
    np.testing.assert_array_equal(moved.faces, mesh.faces)
    np.testing.assert_allclose(
        surface_vertex_normals(moved), surface_vertex_normals(mesh), atol=1e-12
    )

    reflection = np.diag([-1.0, 1.0, 1.0, 1.0])
    reflected = transform_surface_ras(mesh, reflection)
    np.testing.assert_array_equal(reflected.faces, mesh.faces[:, ::-1])
    assert np.mean(
        np.einsum(
            "ij,ij->i",
            surface_vertex_normals(reflected),
            reflected.vertices - np.mean(reflected.vertices, axis=0),
        )
    ) > 0


def test_vtk_round_trip_preserves_thickness_point_data(tmp_path: Path) -> None:
    result = closest_point_thickness(_load_surface("inner"), _load_surface("outer"))
    path = write_legacy_vtk(result.inner_surface, tmp_path / "inner_with_thickness.vtk")

    restored = read_legacy_vtk(path)
    np.testing.assert_allclose(restored.vertices, result.inner_surface.vertices, atol=1e-8)
    np.testing.assert_array_equal(restored.faces, result.inner_surface.faces)
    for name, values in result.inner_surface.point_data.items():
        np.testing.assert_allclose(restored.point_data[name], values, atol=1e-8)


def test_thickness_rejects_surfaces_with_incompatible_data_and_empty_inputs() -> None:
    inner = _load_surface("inner")
    outer = _load_surface("outer")
    with pytest.raises(ValueError, match="one value row"):
        SurfaceMesh(inner.vertices, inner.faces, point_data={"bad": np.ones(2)})

    empty = nib.Nifti1Image(np.ones((4, 4, 4), dtype=np.float32), np.eye(4))
    with pytest.raises(ValueError, match="does not cross"):
        extract_levelset_surface(empty)
