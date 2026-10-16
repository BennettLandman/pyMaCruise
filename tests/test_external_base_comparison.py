"""Exercise portable output mapping and rejection of invalid comparisons."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np
import pytest
from scipy.io import savemat

from pymacruise.surfaces import SurfaceMesh, write_legacy_vtk

SCRIPT = Path(__file__).parents[1] / 'scripts/compare_external_base_run.py'
spec = importlib.util.spec_from_file_location('external_base_comparison', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def case(tmp_path, *, oblique=False, surface_offset=0.0):
    external = tmp_path / 'external/MaCRUISE/Output'
    run = tmp_path / 'python'
    external.mkdir(parents=True)
    run.mkdir()
    affine = np.eye(4)
    if oblique:
        angle = np.deg2rad(8)
        affine[:2, :2] = [[np.cos(angle), -np.sin(angle)],
                          [np.sin(angle), np.cos(angle)]]
    working = np.eye(4)
    working[:3, 3] = nib.affines.apply_affine(affine, [1.5, 1.5, 1.5]) - 1.5
    t1_values = np.arange(64, dtype=np.float32).reshape(4, 4, 4)
    labels = np.zeros((4, 4, 4), dtype=np.int16)
    labels[1:3, 1:3, 1:3] = 100

    def save(path, values, transform):
        path.parent.mkdir(parents=True, exist_ok=True)
        nib.save(nib.Nifti1Image(values, transform), path)
        return path

    t1 = save(tmp_path / 'input/T1.nii.gz', t1_values, affine)
    maseg = save(tmp_path / 'input/MaSeg.nii.gz', labels, affine)
    save(external / 'Raw/T1.nii.gz', t1_values, affine)
    save(external / 'Raw/Seg.nii.gz', labels, affine)
    manifest = {'status': 'complete', 'options': {'regrid_inputs': False},
                'inputs': {'t1': {'path': str(t1), 'sha256': module.sha256(t1)},
                           'maseg': {'path': str(maseg), 'sha256': module.sha256(maseg)}},
                'stages': {}, 'derived_resources': []}

    def record(stage, path):
        manifest['stages'].setdefault(stage, {'status': 'complete', 'outputs': []})
        manifest['stages'][stage]['outputs'].append({'path': str(path), 'sha256': module.sha256(path)})

    for stage in ('n4', 'toads-first', 'dura', 'toads-second', 'base-tissue-preparation'):
        path = run / f'stages/{stage}/attempt-02/log.txt'
        path.parent.mkdir(parents=True)
        path.write_text('fixture\n')
        record(stage, path)
    # A stale attempt must never be selected by the comparator.
    save(run / 'stages/maseg-refinement/attempt-01/SEG_refine.nii.gz', labels + 1, working)
    record('maseg-refinement', save(run / 'stages/maseg-refinement/attempt-02/SEG_refine.nii.gz', labels, working))
    save(external / 'SegRefine/SEG_refine.nii.gz', labels, affine)
    vertices = np.array([[1, 1, 1], [2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=float)
    faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
    for name in ('inner', 'central', 'outer'):
        save(external / f'Surfaces/target_image_GMimg_{name}Vol.nii.gz', t1_values, affine)
        record('cruise-core', save(run / f'stages/cruise-core/attempt-02/{name}.nii.gz', t1_values, working))
        src_mesh = SurfaceMesh(vertices, faces)
        write_legacy_vtk(src_mesh, external / f'Surfaces/target_image_GMimg_{name}Surf.vtk')
        py_mesh = SurfaceMesh(nib.affines.apply_affine(working, vertices) + [surface_offset, 0, 0], faces)
        path = run / f'stages/surface-thickness/attempt-02/{name}_surface.vtk'
        path.parent.mkdir(parents=True, exist_ok=True)
        write_legacy_vtk(py_mesh, path)
        record('surface-thickness', path)
    save(external / 'Thickness/target_image_GMimg_innerVol_thick.nii.gz', t1_values, affine)
    record('surface-thickness', save(run / 'stages/surface-thickness/attempt-02/thickness_levelset.nii.gz', t1_values, working))
    savemat(external / 'Thickness/target_image_GMimg_innerVol_thickness.mat', {'whole_thickness': np.ones(4)})
    path = run / 'stages/surface-thickness/attempt-02/inner_surface_thickness.vtk'
    write_legacy_vtk(SurfaceMesh(nib.affines.apply_affine(working, vertices), faces,
                               point_data={'mean_thickness_mm': np.ones(4)}), path)
    record('surface-thickness', path)
    if oblique:
        for role, values in (('t1', t1_values), ('maseg', labels)):
            path = save(run / f'resources/{role}.nii.gz', values, working)
            manifest['derived_resources'].append({'role': role, 'path': str(path)})
    (run / 'run.json').write_text(json.dumps(manifest))
    info = tmp_path / 'reference.json'
    info.write_text(json.dumps({'image': 'fixture@sha256:test', 'exit_code': 0}))
    log = tmp_path / 'reference.log'
    log.write_text('synthetic comparator fixture, not a real external run\n')
    return SimpleNamespace(external_output=tmp_path / 'external', python_run=run,
                           t1=t1, maseg=maseg, reference_info=info, reference_log=log,
                           source_mesh_frame='voxel', source_volume_grid='raw-t1',
                           report=tmp_path / 'report.json')


def test_identical_products_follow_active_attempt_and_keep_acceptance_unverified(tmp_path):
    result = module.compare(case(tmp_path))
    assert result['segmentation_and_surface_gates_passed']
    assert result['comparisons']['segmentation']['minimum_dice'] == 1.0
    assert 'attempt-02' in result['comparisons']['segmentation']['python_file']
    assert result['overall_port_acceptance'] == 'unverified'


def test_known_surface_offset_fails_mean_distance_gate(tmp_path):
    result = module.compare(case(tmp_path, surface_offset=0.05))
    assert not result['segmentation_and_surface_gates_passed']
    assert result['comparisons']['inner_surface']['mean_absolute_distance_mm'] > 0.01
    assert result['comparisons']['inner_surface']['hausdorff_mm'] < 1


def test_oblique_axis_only_mapping_compares_in_original_physical_frame(tmp_path):
    result = module.compare(case(tmp_path, oblique=True))
    assert result['comparisons']['inner_surface']['hausdorff_mm'] < 1e-6


def test_physical_mode_rejects_oblique_grid_mismatch(tmp_path):
    args = case(tmp_path, oblique=True)
    args.source_volume_grid = 'physical'
    with pytest.raises(ValueError, match='physical grids differ'):
        module.compare(args)


def test_wrong_raw_input_is_rejected(tmp_path):
    args = case(tmp_path)
    path = args.external_output / 'MaCRUISE/Output/Raw/T1.nii.gz'
    nib.save(nib.Nifti1Image(np.zeros((4, 4, 4)), np.eye(4)), path)
    with pytest.raises(ValueError, match='does not preserve'):
        module.compare(args)


def test_mutated_python_output_is_rejected(tmp_path):
    args = case(tmp_path)
    path = args.python_run / 'stages/cruise-core/attempt-02/inner.nii.gz'
    path.write_bytes(b'corrupted')
    with pytest.raises(ValueError, match='hash mismatch'):
        module.compare(args)


def test_missing_external_label_has_zero_dice(tmp_path):
    args = case(tmp_path)
    path = args.external_output / 'MaCRUISE/Output/SegRefine/SEG_refine.nii.gz'
    nib.save(nib.Nifti1Image(np.zeros((4, 4, 4), dtype=np.int16), np.eye(4)), path)
    result = module.compare(args)
    assert result['comparisons']['segmentation']['per_label_dice']['100'] == 0
    assert not result['segmentation_and_surface_gates_passed']


def test_failed_reference_run_is_rejected(tmp_path):
    args = case(tmp_path)
    args.reference_info.write_text(json.dumps({'image': 'test', 'exit_code': 1}))
    with pytest.raises(ValueError, match='successful external exit code'):
        module.compare(args)


def test_cli_writes_error_report_for_missing_external_products(tmp_path):
    args = case(tmp_path)
    (args.external_output / 'MaCRUISE/Output/SegRefine/SEG_refine.nii.gz').unlink()
    argv = []
    for key, value in vars(args).items():
        argv += ['--' + key.replace('_', '-'), str(value)]
    assert module.main(argv) == 2
    assert 'comparison_error' in json.loads(args.report.read_text())
