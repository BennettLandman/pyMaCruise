from pathlib import Path
from contextlib import nullcontext
import json

import nibabel as nib
import numpy as np
import pytest

from pymacruise.base_native_export import (
    NATIVE_OUTPUT_NAMES, export_base_native_outputs, restore_image_to_original_grid,
)
from pymacruise.image_io import orient_to_ras_without_obliquity, regrid_to_orthogonal_ras
from pymacruise.surfaces import SurfaceMesh, read_legacy_vtk, write_legacy_vtk


def original_image(*, nonras=False):
    angle = np.deg2rad(13)
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    axes = np.array([[0, -1, 0], [1, 0, 0], [0, 0, -1]]) if nonras else np.eye(3)
    affine = np.eye(4)
    affine[:3, :3] = rotation @ axes @ np.diag([.8, .9, 1.2])
    affine[:3, 3] = [-4, 2, 10]
    data = np.arange(7 * 8 * 9, dtype=np.uint16).reshape(7, 8, 9)
    image = nib.Nifti1Image(data, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=2)
    return image


def assert_native(image, original):
    assert image.shape == original.shape
    np.testing.assert_allclose(image.affine, original.affine, rtol=0, atol=1e-6)
    np.testing.assert_allclose(image.get_qform(), original.get_qform(), rtol=0, atol=1e-6)
    np.testing.assert_allclose(image.get_sform(), original.get_sform(), rtol=0, atol=1e-6)
    assert int(image.header['qform_code']) == int(original.header['qform_code'])
    assert int(image.header['sform_code']) == int(original.header['sform_code'])


def test_oblique_native_export_samples_physical_linear_field_and_discrete_labels(tmp_path):
    original = original_image()
    working = regrid_to_orthogonal_ras(original)
    points = np.moveaxis(np.indices(working.shape), 0, -1)
    world = nib.affines.apply_affine(working.affine, points)
    plane = 2 * world[..., 0] + world[..., 1] - .5 * world[..., 2]
    continuous = nib.Nifti1Image(plane.astype(np.float32), working.affine)
    restored = restore_image_to_original_grid(continuous, original, working)
    native_points = np.moveaxis(np.indices(original.shape), 0, -1)
    native_world = nib.affines.apply_affine(original.affine, native_points)
    expected = 2 * native_world[..., 0] + native_world[..., 1] - .5 * native_world[..., 2]
    np.testing.assert_allclose(np.asarray(restored.dataobj), expected, rtol=0, atol=1e-5)
    labels = np.where(world[..., 0] > -2, 100, 45).astype(np.int32)
    segmentation = restore_image_to_original_grid(
        nib.Nifti1Image(labels, working.affine), original, working, labels=True,
    )
    assert set(np.unique(np.asarray(segmentation.dataobj))) == {45, 100}
    for name, image in [('continuous', restored), ('labels', segmentation)]:
        path = tmp_path / f'{name}.nii.gz'
        nib.save(image, path)
        saved = nib.load(path)
        assert_native(saved, original)
    assert np.issubdtype(nib.load(tmp_path / 'labels.nii.gz').get_data_dtype(), np.integer)


def test_orientation_only_export_reverses_permutations_and_flips_without_interpolation():
    original = original_image(nonras=True)
    with pytest.warns(UserWarning):
        working = orient_to_ras_without_obliquity(original)
    restored = restore_image_to_original_grid(
        working, original, working, orientation_only=True, labels=True,
    )
    np.testing.assert_array_equal(np.asarray(restored.dataobj), np.asarray(original.dataobj))
    assert_native(restored, original)


@pytest.mark.parametrize('orientation_only', [False, True])
def test_complete_native_export_retains_grid_surface_geometry_and_thickness(tmp_path, orientation_only):
    original = original_image(nonras=orientation_only)
    if orientation_only:
        with pytest.warns(UserWarning):
            working = orient_to_ras_without_obliquity(original)
    else:
        working = regrid_to_orthogonal_ras(original)
    def save(name, image):
        path = tmp_path / name
        nib.save(image, path)
        return path
    reference_path = save('reference.nii.gz', original)
    working_path = save('working.nii.gz', working)
    levelset = nib.Nifti1Image(np.full(working.shape, -1, dtype=np.float32), working.affine)
    levelset_path = save('levelset.nii.gz', levelset)
    source_dir = tmp_path / 'surfaces'
    source_dir.mkdir()
    nib.save(levelset, source_dir / 'thickness_levelset.nii.gz')
    (source_dir / 'thickness_summary.json').write_text(json.dumps({'mean_mm': 2.5}))
    voxel_vertices = np.array([[2, 2, 2], [3, 2, 2], [2, 3, 2], [2, 2, 3]], dtype=float)
    vertices = nib.affines.apply_affine(working.affine, voxel_vertices)
    faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
    mesh = SurfaceMesh(vertices, faces, point_data={'mean_thickness_mm': np.arange(4, dtype=float)})
    for name in ('inner_surface.vtk', 'central_surface.vtk', 'outer_surface.vtk', 'inner_surface_thickness.vtk'):
        write_legacy_vtk(mesh, source_dir / name)
    with pytest.warns(UserWarning) if orientation_only else nullcontext():
        paths = export_base_native_outputs(
            reference_path=reference_path, working_reference_path=working_path,
            segmentation_path=working_path, inner_path=levelset_path, central_path=levelset_path,
            outer_path=levelset_path, surface_dir=source_dir, output_dir=tmp_path / 'native',
            orientation_only=orientation_only,
        )
    assert set(paths) == set(NATIVE_OUTPUT_NAMES)
    for name, path in paths.items():
        if name.endswith('.nii.gz'):
            assert_native(nib.load(path), original)
    assert np.issubdtype(nib.load(paths['SEG_refine.nii.gz']).get_data_dtype(), np.integer)
    final_mesh = read_legacy_vtk(paths['inner_surface_thickness.vtk'])
    expected_affine = nib.as_closest_canonical(original).affine if orientation_only else working.affine
    np.testing.assert_allclose(final_mesh.vertices, nib.affines.apply_affine(expected_affine, voxel_vertices), atol=1e-6)
    np.testing.assert_array_equal(final_mesh.faces, faces)
    np.testing.assert_array_equal(final_mesh.point_data['mean_thickness_mm'], np.arange(4))


def test_native_export_rejects_fractional_labels_and_wrong_working_grid():
    original = original_image()
    fractional = nib.Nifti1Image(np.full(original.shape, .5), original.affine)
    with pytest.raises(ValueError, match='integer labels'):
        restore_image_to_original_grid(fractional, original, original, labels=True)
    with pytest.raises(ValueError, match='working T1 grid'):
        restore_image_to_original_grid(fractional, original, regrid_to_orthogonal_ras(original))
