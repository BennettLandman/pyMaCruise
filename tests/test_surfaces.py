from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.surfaces import (
    SurfaceMesh,
    convert_surface_to_ras,
    extract_levelset_surface,
    read_legacy_vtk,
    write_legacy_vtk,
)

SURFACE_FIXTURES = Path(__file__).parent / "data" / "synthetic" / "surfaces"


def _mesh_topology(mesh: SurfaceMesh) -> tuple[int, int, int]:
    edges = np.sort(
        np.concatenate(
            (mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]])
        ),
        axis=1,
    )
    unique_edges, edge_counts = np.unique(edges, axis=0, return_counts=True)
    assert np.all(edge_counts == 2)
    return len(mesh.vertices) - len(unique_edges) + len(mesh.faces), len(unique_edges), len(mesh.faces)


@pytest.mark.parametrize(
    ("name", "radii"),
    [("sphere", (5.5, 5.5, 5.5)), ("ellipsoid", (5.5, 4.5, 6.0))],
)
def test_synthetic_sphere_and_ellipsoid_extract_closed_ras_mesh(name, radii) -> None:
    image = nib.load(SURFACE_FIXTURES / f"{name}.nii.gz")
    mesh = extract_levelset_surface(image)

    assert mesh.coordinate_frame == "RAS"
    assert mesh.units == "mm"
    assert np.isfinite(mesh.vertices).all()
    assert _mesh_topology(mesh)[0] == 2
    assert np.isclose(
        np.sum(
            np.einsum(
                "ij,ij->i",
                mesh.vertices[mesh.faces][:, 0],
                np.cross(mesh.vertices[mesh.faces][:, 1], mesh.vertices[mesh.faces][:, 2]),
            )
        ),
        6 * (4 * np.pi / 3) * np.prod(radii),
        rtol=0.08,
    )
    center = image.affine[:3, :3] @ ((np.asarray(image.shape) - 1) / 2)
    normalized_radius = np.sqrt(
        np.sum(((mesh.vertices - center) / np.asarray(radii)) ** 2, axis=1)
    )
    assert np.max(np.abs(normalized_radius - 1.0)) < 0.04


def test_oblique_levelset_affine_is_applied_to_mesh_coordinates() -> None:
    image = nib.load(SURFACE_FIXTURES / "sphere.nii.gz")
    affine = image.affine.copy()
    angle = np.deg2rad(23)
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]]
    )
    affine[:3, :3] = rotation @ affine[:3, :3]
    affine[:3, 3] = [12, -7, 4]
    oblique = nib.Nifti1Image(np.asarray(image.dataobj), affine)

    mesh = extract_levelset_surface(oblique)
    center = affine[:3, :3] @ ((np.asarray(image.shape) - 1) / 2) + affine[:3, 3]
    distances = np.linalg.norm(mesh.vertices - center, axis=1)
    assert np.max(np.abs(distances - 5.5)) < 0.12
    assert np.allclose(mesh.vertices.mean(axis=0), center, atol=0.2)


def test_non_ras_levelset_warns_and_generates_same_ras_surface() -> None:
    ras = nib.load(SURFACE_FIXTURES / "sphere.nii.gz")
    ras_mesh = extract_levelset_surface(ras)
    lps_affine = ras.affine.copy()
    lps_affine[:2, :2] *= -1
    lps_affine[0, 3] = (ras.shape[0] - 1) * 0.5
    lps_affine[1, 3] = (ras.shape[1] - 1) * 0.5
    lps_image = nib.Nifti1Image(np.asarray(ras.dataobj)[::-1, ::-1], lps_affine)

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        lps_mesh = extract_levelset_surface(lps_image)
    np.testing.assert_allclose(
        np.sort(ras_mesh.vertices, axis=0), np.sort(lps_mesh.vertices, axis=0), atol=1e-6
    )


def test_lps_surface_conversion_warns_flips_coordinates_and_winding() -> None:
    mesh = extract_levelset_surface(nib.load(SURFACE_FIXTURES / "sphere.nii.gz"))
    with pytest.warns(UserWarning, match="LPS.*converting"):
        converted = convert_surface_to_ras(mesh, input_frame="LPS")
    with pytest.warns(UserWarning, match="LPS.*converting"):
        converted_again = convert_surface_to_ras(converted, input_frame="LPS")
    np.testing.assert_allclose(converted.vertices, -mesh.vertices * [1, 1, -1])
    np.testing.assert_array_equal(converted.faces, mesh.faces[:, ::-1])
    np.testing.assert_allclose(converted_again.vertices, mesh.vertices)


def test_vtk_ascii_surface_roundtrip_and_legacy_lps_frame(tmp_path: Path) -> None:
    mesh = extract_levelset_surface(nib.load(SURFACE_FIXTURES / "sphere.nii.gz"))
    output = write_legacy_vtk(mesh, tmp_path / "sphere.vtk")
    restored = read_legacy_vtk(output)
    np.testing.assert_array_equal(restored.vertices, mesh.vertices)
    np.testing.assert_array_equal(restored.faces, mesh.faces)

    lps_path = tmp_path / "lps.vtk"
    lps_path.write_text(
        "# vtk DataFile Version 3.0\nOld mesh\nASCII\nDATASET POLYDATA\n"
        "POINTS 4 float\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n"
        "POLYGONS 4 16\n3 0 2 1\n3 0 1 3\n3 0 3 2\n3 1 2 3\n",
        encoding="ascii",
    )
    with pytest.raises(ValueError, match="requires input_frame"):
        read_legacy_vtk(lps_path)
    with pytest.warns(UserWarning, match="LPS.*converting"):
        converted = read_legacy_vtk(lps_path, input_frame="LPS")
    np.testing.assert_array_equal(converted.vertices[1], [-1, 0, 0])
    assert _mesh_topology(converted)[0] == 2


def test_surface_extraction_rejects_non_crossing_and_nonfinite_images() -> None:
    affine = np.eye(4)
    no_crossing = nib.Nifti1Image(np.ones((5, 5, 5), dtype=np.float32), affine)
    with pytest.raises(ValueError, match="does not cross"):
        extract_levelset_surface(no_crossing)
    invalid = np.ones((5, 5, 5), dtype=np.float32)
    invalid[2, 2, 2] = np.nan
    with pytest.raises(ValueError, match="NaN or infinite"):
        extract_levelset_surface(nib.Nifti1Image(invalid, affine))
