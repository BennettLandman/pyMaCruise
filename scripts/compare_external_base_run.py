#!/usr/bin/env python3
"""Compare independent external BASE products with a restartable Python run.

The external layout is the native Output tree documented by MASILab's v3.1
container. Coordinate contracts are explicit rather than inferred from mesh
bounds. This computes product agreement, not a whole-port acceptance claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
from scipy.io import loadmat

from pymacruise.surfaces import read_legacy_vtk
from pymacruise.validation import (
    dice_coefficient, label_dice_scores, mesh_surface_distance_metrics,
    voxel_error_metrics,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def image(path: Path):
    if not path.is_file() or path.is_symlink():
        raise ValueError(f'Missing regular image: {path}')
    result = nib.load(path)
    values = np.asanyarray(result.dataobj)
    if (len(result.shape) != 3 or not np.issubdtype(values.dtype, np.number)
            or np.iscomplexobj(values) or not np.isfinite(values).all()
            or not np.isfinite(result.affine).all()):
        raise ValueError(f'Expected a finite 3D image: {path}')
    return result, values


def summary(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=np.float64).ravel()
    if not len(values) or not np.isfinite(values).all():
        raise ValueError('Thickness samples must be nonempty and finite')
    return {'count': len(values), 'mean_mm': float(values.mean()),
            'std_mm': float(values.std()), 'min_mm': float(values.min()),
            'median_mm': float(np.median(values)), 'max_mm': float(values.max())}


def find_output(root: Path) -> Path:
    candidates = {p.parent.parent for p in root.rglob('SEG_refine.nii.gz')
                  if p.parent.name == 'SegRefine'
                  and (p.parent.parent / 'Raw/T1.nii.gz').is_file()
                  and (p.parent.parent / 'Surfaces').is_dir()}
    if len(candidates) != 1:
        raise ValueError(f'Expected one native external Output tree; found {len(candidates)} under {root}')
    return candidates.pop()


def compare(args) -> dict:
    source = find_output(args.external_output.resolve())
    run_path = args.python_run.resolve() / 'run.json'
    manifest = json.loads(run_path.read_text())
    required_stages = {'n4', 'toads-first', 'dura', 'toads-second',
                       'base-tissue-preparation', 'cruise-core',
                       'maseg-refinement', 'surface-thickness'}
    if manifest.get('status') != 'complete' or set(manifest['stages']) != required_stages:
        raise ValueError('Python manifest must contain all eight completed BASE stages')
    verified = 0
    for name, stage in manifest['stages'].items():
        if stage.get('status') != 'complete':
            raise ValueError(f'Incomplete Python stage: {name}')
        for record in stage['outputs']:
            path = Path(record['path'])
            if not path.is_file() or path.is_symlink() or sha256(path) != record['sha256']:
                raise ValueError(f'Python output hash mismatch or missing file: {path}')
            verified += 1

    def active(stage: str, filename: str) -> Path:
        paths = [Path(r['path']) for r in manifest['stages'][stage]['outputs']
                 if Path(r['path']).name == filename]
        if len(paths) != 1:
            raise ValueError(f'Expected one active {stage}/{filename}')
        return paths[0]

    info = json.loads(args.reference_info.read_text())
    if info.get('exit_code') != 0 or not info.get('image'):
        raise ValueError('Reference record must identify the image and a successful external exit code')
    if not args.reference_log.is_file() or args.reference_log.stat().st_size == 0:
        raise ValueError('A nonempty external launch log is required')

    raw_t1, raw_values = image(source / 'Raw/T1.nii.gz')
    raw_seg, raw_labels = image(source / 'Raw/Seg.nii.gz')
    t1, t1_values = image(args.t1)
    maseg, maseg_values = image(args.maseg)
    for name, a, b, av, bv in (
        ('T1', raw_t1, t1, raw_values, t1_values),
        ('MaSeg', raw_seg, maseg, raw_labels, maseg_values),
    ):
        if a.shape != b.shape or not np.array_equal(av, bv) or not np.allclose(a.affine, b.affine, atol=1e-4, rtol=0):
            raise ValueError(f'External Raw/{name} does not preserve the declared input grid and values')
    if raw_t1.shape != raw_seg.shape or not np.allclose(raw_t1.affine, raw_seg.affine, atol=1e-4, rtol=0):
        raise ValueError('External raw T1 and MaSeg grids differ')
    if nib.aff2axcodes(t1.affine) != ('R', 'A', 'S'):
        raise ValueError('This native-output adapter requires the supplied RAS input packet')

    # Prove the Python run used these inputs, not another same-shaped scan.
    for name, supplied, supplied_path in (('T1', t1_values, args.t1),
                                          ('MaSeg', maseg_values, args.maseg)):
        resource = next((r for r in manifest.get('derived_resources', [])
                         if r.get('role') == ('t1' if name == 'T1' else 'maseg')), None)
        original = manifest['inputs'].get('t1' if name == 'T1' else 'maseg')
        if not isinstance(original, dict) or original.get('sha256') != sha256(supplied_path):
            raise ValueError(f'Python original {name} hash differs from the supplied input')
        if resource is not None:
            candidate_path = Path(resource['path'])
        elif isinstance(original, dict):
            candidate_path = Path(original['path'])
        elif isinstance(original, str):
            candidate_path = Path(original)
        else:
            raise ValueError(f'Missing Python input provenance for {name}')
        _, candidate_values = image(candidate_path)
        if not np.array_equal(supplied, candidate_values):
            raise ValueError(f'Python {name} values differ from the supplied input')
    if manifest.get('options', {}).get('regrid_inputs', False):
        raise ValueError('Use the documented no-regrid Python recipe for this adapter')

    report = {
        'scope': 'independent BASE product comparison; preprocessing contracts may differ',
        'overall_port_acceptance': 'unverified',
        'reference': info,
        'reference_info_sha256': sha256(args.reference_info),
        'reference_log_sha256': sha256(args.reference_log),
        'external_output': str(source), 'python_run': str(run_path.parent),
        'python_manifest_sha256': sha256(run_path),
        'verified_python_stage_records': verified,
        'source_mesh_frame': args.source_mesh_frame,
        'source_volume_grid': args.source_volume_grid,
        'raw_t1_affine_ras_mm': raw_t1.affine.tolist(),
        'input_identity': {'T1': sha256(args.t1), 'MaSeg': sha256(args.maseg)},
        'comparisons': {}, 'limitations': [
            'External skull stripping and Python supplied-mask preprocessing differ unless independently shown equivalent.',
            'Scalar and thickness acceptance limits have not been established; metrics do not imply an overall pass.',
            'Raw-grid mode is an explicit source export contract; it does not independently prove the output voxel ordering.',
        ],
    }
    pairs = [('segmentation', source / 'SegRefine/SEG_refine.nii.gz',
              active('maseg-refinement', 'SEG_refine.nii.gz'))]
    for name in ('inner', 'central', 'outer'):
        pairs.append((name, source / f'Surfaces/target_image_GMimg_{name}Vol.nii.gz',
                      active('cruise-core', f'{name}.nii.gz')))
    pairs.append(('thickness_volume', source / 'Thickness/target_image_GMimg_innerVol_thick.nii.gz',
                  active('surface-thickness', 'thickness_levelset.nii.gz')))
    python_affine = None
    for name, src_path, py_path in pairs:
        src_img, src_values = image(src_path)
        py_img, py_values = image(py_path)
        if src_img.shape != raw_t1.shape or py_img.shape != raw_t1.shape:
            raise ValueError(f'{name} is not on the full native input grid')
        if nib.aff2axcodes(py_img.affine) != ('R', 'A', 'S'):
            raise ValueError(f'Python {name} is not RAS')
        if args.source_volume_grid == 'physical':
            if (not np.allclose(src_img.affine, raw_t1.affine, atol=1e-4, rtol=0)
                    or not np.allclose(src_img.affine, py_img.affine, atol=1e-4, rtol=0)):
                raise ValueError(f'{name} physical grids differ; do not compare arrays by shape alone')
        else:
            expected_affine = np.diag([*nib.affines.voxel_sizes(raw_t1.affine), 1.0])
            center = (np.asarray(raw_t1.shape) - 1) / 2
            expected_affine[:3, 3] = (nib.affines.apply_affine(raw_t1.affine, center)
                                       - expected_affine[:3, :3] @ center)
            if not (np.allclose(py_img.affine, expected_affine, atol=1e-4, rtol=0)
                    or np.allclose(py_img.affine, raw_t1.affine, atol=1e-4, rtol=0)):
                raise ValueError(f'Python {name} does not have the declared axis-only native mapping')
        entry = {'external_file': str(src_path), 'external_sha256': sha256(src_path),
                 'python_file': str(py_path), 'python_sha256': sha256(py_path),
                 'external_affine': src_img.affine.tolist(), 'python_affine': py_img.affine.tolist()}
        if name == 'segmentation':
            if not np.array_equal(src_values, np.rint(src_values)) or not np.array_equal(py_values, np.rint(py_values)):
                raise ValueError('Segmentation label images must contain integer values')
            scores = label_dice_scores(src_values, py_values)
            if not scores:
                raise ValueError('Refined segmentations contain no nonbackground labels')
            entry.update(per_label_dice={str(k): v for k, v in scores.items()},
                         per_label_voxel_counts={str(k): {
                             'external': int(np.count_nonzero(src_values == k)),
                             'python': int(np.count_nonzero(py_values == k))} for k in scores},
                         differing_voxels=int(np.count_nonzero(src_values != py_values)),
                         minimum_dice=min(scores.values()),
                         labels_below_threshold=[int(k) for k, v in scores.items() if v < 0.97],
                         passed=all(v >= 0.97 for v in scores.values()))
        else:
            entry['voxel_error'] = voxel_error_metrics(src_values, py_values)
            if name in ('inner', 'central', 'outer'):
                entry['nonpositive_mask_dice'] = dice_coefficient(src_values <= 0, py_values <= 0)
                python_affine = py_img.affine
            else:
                union = (src_values > 0) | (py_values > 0)
                entry['positive_mask_dice'] = dice_coefficient(src_values > 0, py_values > 0)
                entry['positive_union_error'] = (voxel_error_metrics(src_values[union], py_values[union])
                                                 if union.any() else None)
        report['comparisons'][name] = entry

    for name in ('inner', 'central', 'outer'):
        src_path = source / f'Surfaces/target_image_GMimg_{name}Surf.vtk'
        py_path = active('surface-thickness', f'{name}_surface.vtk')
        src_mesh = read_legacy_vtk(src_path, input_frame='RAS')
        py_mesh = read_legacy_vtk(py_path, input_frame='RAS')
        src_vertices = src_mesh.vertices
        py_vertices = py_mesh.vertices
        if args.source_mesh_frame == 'voxel':
            src_vertices = nib.affines.apply_affine(raw_t1.affine, src_vertices)
        if args.source_volume_grid == 'raw-t1':
            py_vertices = nib.affines.apply_affine(raw_t1.affine @ np.linalg.inv(python_affine), py_vertices)
        print(f'Comparing {name} triangle surfaces...', flush=True)
        distances = mesh_surface_distance_metrics(src_vertices, src_mesh.faces, py_vertices, py_mesh.faces)
        report['comparisons'][f'{name}_surface'] = {
            'external_file': str(src_path), 'external_sha256': sha256(src_path),
            'python_file': str(py_path), 'python_sha256': sha256(py_path),
            **distances,
            'passed': distances['mean_absolute_distance_mm'] <= 0.01 and distances['hausdorff_mm'] <= 1.0,
        }
    mat_path = source / 'Thickness/target_image_GMimg_innerVol_thickness.mat'
    if not mat_path.is_file():
        raise ValueError(f'Missing external surface thickness MAT file: {mat_path}')
    scalar_mesh = read_legacy_vtk(active('surface-thickness', 'inner_surface_thickness.vtk'), input_frame='RAS')
    src_sample = loadmat(mat_path)['whole_thickness'].ravel()
    py_sample = scalar_mesh.point_data['mean_thickness_mm']
    src_inner = read_legacy_vtk(source / 'Surfaces/target_image_GMimg_innerSurf.vtk', input_frame='RAS')
    if len(src_sample) != len(src_inner.vertices):
        raise ValueError('External surface thickness sample count differs from the inner mesh')
    report['comparisons']['surface_thickness'] = {
        'external_file': str(mat_path), 'external_sha256': sha256(mat_path),
        'external_distribution': summary(src_sample), 'python_distribution': summary(py_sample),
        'vertexwise_comparison': 'unverified; independently generated meshes have no established vertex correspondence',
    }
    report['segmentation_and_surface_gates_passed'] = all(
        report['comparisons'][name]['passed']
        for name in ('segmentation', 'inner_surface', 'central_surface', 'outer_surface'))
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('external-output', 'python-run', 't1', 'maseg', 'reference-info', 'reference-log', 'report'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--source-mesh-frame', choices=('voxel', 'ras-mm'), required=True)
    parser.add_argument('--source-volume-grid', choices=('physical', 'raw-t1'), default='physical')
    args = parser.parse_args(argv)
    try:
        report = compare(args)
        status = 0 if report['segmentation_and_surface_gates_passed'] else 1
    except (ValueError, KeyError, OSError) as exc:
        report = {'overall_port_acceptance': 'unverified', 'comparison_error': str(exc)}
        status = 2
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(f'Report: {args.report}; exit code: {status}', flush=True)
    return status


if __name__ == '__main__':
    sys.exit(main())
