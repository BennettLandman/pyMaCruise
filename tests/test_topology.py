from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest
from scipy import ndimage

from pymacruise.topology import topology_correct, topology_numbers_618


def test_618_topology_numbers_cover_regular_and_critical_neighborhoods() -> None:
    empty = np.zeros((3, 3, 3), dtype=bool)
    assert topology_numbers_618(empty) == (0, 1)

    one_face = empty.copy()
    one_face[0, 1, 1] = True
    assert topology_numbers_618(one_face) == (1, 1)

    two_disconnected_faces = empty.copy()
    two_disconnected_faces[0, 1, 1] = True
    two_disconnected_faces[2, 1, 1] = True
    assert topology_numbers_618(two_disconnected_faces) == (2, 1)

    surrounded = np.ones((3, 3, 3), dtype=bool)
    surrounded[1, 1, 1] = False
    assert topology_numbers_618(surrounded) == (1, 0)


def test_618_foreground_ignores_components_not_touching_a_face_neighbor() -> None:
    neighborhood = np.zeros((3, 3, 3), dtype=bool)
    neighborhood[0, 1, 0] = True  # isolated edge-adjacent component

    assert topology_numbers_618(neighborhood) == (0, 1)


def test_fast_marching_keeps_separate_painted_objects_separate_and_is_repeatable() -> None:
    shape = (15, 15, 15)
    affine = np.diag([0.9, 1.1, 1.7, 1.0])
    affine[:3, 3] = [-12.0, 8.0, -22.0]
    values = np.full(shape, 255, dtype=np.float32)
    values[0, 0, 0] = 0
    paint = np.zeros(shape, dtype=np.uint8)
    paint[5, 7, 7] = 1
    paint[9, 7, 7] = 1
    image = nib.Nifti1Image(values, affine)
    paint_image = nib.Nifti1Image(paint, affine)

    corrected = topology_correct(image, paint_image).corrected
    repeated = topology_correct(image, paint_image).corrected
    corrected_data = np.asarray(corrected.dataobj)
    assert corrected_data.dtype == np.float32
    np.testing.assert_array_equal(corrected_data, np.asarray(repeated.dataobj))
    assert corrected.shape == shape
    assert nib.aff2axcodes(corrected.affine) == ("R", "A", "S")
    np.testing.assert_allclose(corrected.affine, affine)

    _, components = ndimage.label(
        corrected_data > 127.5,
        structure=ndimage.generate_binary_structure(3, 1),
    )
    assert components == 2


@pytest.mark.parametrize("shape_kind, expected_euler", [("torus", 0), ("shell", 2)])
def test_seeded_handles_and_cavity_keep_their_cubical_euler_characteristic(
    shape_kind: str,
    expected_euler: int,
) -> None:
    size = 17
    x, y, z = np.indices((size, size, size))
    radius = np.sqrt((x - 8) ** 2 + (y - 8) ** 2 + (z - 8) ** 2)
    if shape_kind == "torus":
        radial = np.sqrt((x - 8) ** 2 + (y - 8) ** 2)
        object_mask = (radial - 4) ** 2 + (z - 8) ** 2 <= 1.8**2
    else:
        object_mask = (radius >= 4) & (radius <= 6)
    intensity = object_mask.astype(np.float32) * 255
    # The zero corner makes the padded-source intensity range non-degenerate.
    intensity[0, 0, 0] = 0
    result = topology_correct(
        nib.Nifti1Image(intensity, np.eye(4)),
        nib.Nifti1Image(object_mask.astype(np.uint8), np.eye(4)),
    )
    corrected_mask = np.asarray(result.corrected.dataobj) > 127.5

    assert _cubical_euler(object_mask) == expected_euler
    assert _cubical_euler(corrected_mask) == expected_euler
    np.testing.assert_array_equal(corrected_mask, object_mask)


def test_v3_atlas_refinement_is_a_separate_output_and_preserves_geometry() -> None:
    shape = (9, 9, 9)
    affine = np.eye(4)
    values = np.full(shape, 255, dtype=np.float32)
    values[0, 0, 0] = 0
    values[3, 3, 3] = 0
    paint = np.zeros(shape, dtype=np.uint8)
    paint[4, 4, 4] = 1
    atlas_data = np.zeros(shape, dtype=np.uint8)
    atlas_data[2:5, 2:5, 2:5] = 1

    result = topology_correct(
        nib.Nifti1Image(values, affine),
        nib.Nifti1Image(paint, affine),
        multi_atlas_wm=nib.Nifti1Image(atlas_data, affine),
    )

    assert result.atlas_refined is not None
    assert np.asarray(result.corrected.dataobj)[3, 3, 3] < 200
    assert np.asarray(result.atlas_refined.dataobj)[3, 3, 3] == 255
    np.testing.assert_array_equal(result.atlas_refined.affine, affine)
    assert nib.aff2axcodes(result.atlas_refined.affine) == ("R", "A", "S")


def test_topology_correct_warns_on_non_ras_inputs_and_checks_grids() -> None:
    shape = (9, 10, 11)
    values = np.ones(shape, dtype=np.float32)
    values[0, 0, 0] = 0
    paint = np.zeros(shape, dtype=np.uint8)
    paint[4, 5, 5] = 1
    non_ras = np.diag([-1.0, -1.0, 1.0, 1.0])

    with pytest.warns(UserWarning, match="reorienting the image to RAS") as caught:
        result = topology_correct(
            nib.Nifti1Image(values, non_ras),
            nib.Nifti1Image(paint, non_ras),
        )
    assert len(caught) == 2
    assert nib.aff2axcodes(result.corrected.affine) == ("R", "A", "S")

    shifted = np.eye(4)
    shifted[0, 3] = 1
    with pytest.raises(ValueError, match="same physical grid"):
        topology_correct(
            nib.Nifti1Image(values, np.eye(4)),
            nib.Nifti1Image(paint, shifted),
        )


def test_topology_correct_rejects_empty_paint_and_nonfinite_values() -> None:
    shape = (5, 5, 5)
    values = np.ones(shape, dtype=np.float32)
    empty = np.zeros(shape, dtype=np.uint8)
    image = nib.Nifti1Image(values, np.eye(4))
    with pytest.raises(ValueError, match="no voxels"):
        topology_correct(image, nib.Nifti1Image(empty, np.eye(4)))

    values[0, 0, 0] = np.nan
    paint = empty.copy()
    paint[2, 2, 2] = 1
    with pytest.raises(ValueError, match="NaN"):
        topology_correct(
            nib.Nifti1Image(values, np.eye(4)),
            nib.Nifti1Image(paint, np.eye(4)),
        )


def _cubical_euler(mask: np.ndarray) -> int:
    """Euler characteristic of the union of closed unit voxel cubes."""
    cubes = np.asarray(mask, dtype=bool)
    nx, ny, nz = cubes.shape
    vertices = np.zeros((nx + 1, ny + 1, nz + 1), dtype=bool)
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                vertices[dx : dx + nx, dy : dy + ny, dz : dz + nz] |= cubes

    faces_x = np.zeros((nx + 1, ny, nz), dtype=bool)
    faces_x[:-1] |= cubes
    faces_x[1:] |= cubes
    faces_y = np.zeros((nx, ny + 1, nz), dtype=bool)
    faces_y[:, :-1] |= cubes
    faces_y[:, 1:] |= cubes
    faces_z = np.zeros((nx, ny, nz + 1), dtype=bool)
    faces_z[:, :, :-1] |= cubes
    faces_z[:, :, 1:] |= cubes

    edges_x = np.zeros((nx, ny + 1, nz + 1), dtype=bool)
    edges_y = np.zeros((nx + 1, ny, nz + 1), dtype=bool)
    edges_z = np.zeros((nx + 1, ny + 1, nz), dtype=bool)
    for dy in (0, 1):
        for dz in (0, 1):
            edges_x[:, dy : dy + ny, dz : dz + nz] |= cubes
    for dx in (0, 1):
        for dz in (0, 1):
            edges_y[dx : dx + nx, :, dz : dz + nz] |= cubes
    for dx in (0, 1):
        for dy in (0, 1):
            edges_z[dx : dx + nx, dy : dy + ny, :] |= cubes
    return int(
        vertices.sum()
        - edges_x.sum()
        - edges_y.sum()
        - edges_z.sum()
        + faces_x.sum()
        + faces_y.sum()
        + faces_z.sum()
        - cubes.sum()
    )
