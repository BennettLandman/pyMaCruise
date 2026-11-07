#!/usr/bin/env python3
"""Verify a terminal Python PLUS tree without requiring an original run."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import nibabel as nib
import numpy as np
from pymacruise.plus_full_workflow import PLUS_STAGES
from pymacruise.surfaces import read_legacy_vtk
from compare_external_plus_run import _python_products, _records, sha256, verify_hash_records


def audit(root: Path) -> dict:
    manifest_path = root.resolve() / 'plus-full-workflow.json'
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('status') != 'complete' or manifest.get('terminal') is not True:
        raise ValueError('PLUS workflow is not complete and terminal')
    expected = list(PLUS_STAGES)
    if not manifest.get('settings', {}).get('atlasem_dura_enabled'):
        expected.remove('atlasem_dura')
    if manifest.get('stage_order') != expected or set(manifest.get('stages', {})) != set(expected):
        raise ValueError('Missing, extra or reordered PLUS stages')
    documents = [manifest]
    for name in expected:
        stage = manifest['stages'][name]
        if stage.get('status') != 'complete':
            raise ValueError(f'Incomplete PLUS stage: {name}')
        path = Path(stage['manifest']['path'])
        child = json.loads(path.read_text())
        if child.get('status') != 'complete':
            raise ValueError(f'Incomplete child manifest: {path}')
        documents.append(child)
    surface_path = Path(manifest['stages']['plus_surface_workflow']['manifest']['path'])
    surface = json.loads(surface_path.read_text())
    for name, child in surface['stages'].items():
        if child.get('status') != 'complete':
            raise ValueError(f'Incomplete surface substage: {name}')
    paths = {}
    for document in documents:
        # AtlasEM payloads also use a path->digest mapping rather than
        # ordinary file records. Check that mapping as well as records.
        verify_hash_records(document, 'Python PLUS audit')
        for record in _records(document):
            path = Path(record['path'])
            digest = record['sha256']
            if path in paths and paths[path] != digest:
                raise ValueError(f'Inconsistent file hashes: {path}')
            if path not in paths:
                if not path.is_file() or path.is_symlink() or sha256(path) != digest:
                    raise ValueError(f'Output/input hash verification failed: {path}')
                paths[path] = digest
    images = []
    nonfinite_diagnostic_images = []
    for path in paths:
        if not path.name.endswith('.nii.gz'):
            continue
        image = nib.load(path)
        if nib.aff2axcodes(image.affine) != ('R', 'A', 'S'):
            raise ValueError(f'Non-RAS image: {path}')
        # Energy/state fields can legitimately preserve source IEEE NaNs.
        # Final scientific products are required finite below.
        finite = all(np.isfinite(np.asarray(image.dataobj[x:x+8])).all()
                     for x in range(0, image.shape[0], 8))
        (images if finite else nonfinite_diagnostic_images).append(str(path))
    products = _python_products(manifest)
    meshes = {}
    segmentation = nib.load(products['segmentation'])
    for name, path in products.items():
        if path.name.endswith('.nii.gz'):
            image = nib.load(path)
            if image.shape != segmentation.shape or not np.allclose(image.affine, segmentation.affine, atol=1e-4, rtol=0):
                raise ValueError(f'Final products do not share one full grid: {name}')
            if str(path) in nonfinite_diagnostic_images:
                raise ValueError(f'Nonfinite final scientific product: {path}')
            if name == 'segmentation':
                for x in range(0, image.shape[0], 8):
                    values = np.asarray(image.dataobj[x:x+8])
                    if not np.array_equal(values, np.rint(values)):
                        raise ValueError('Final segmentation is not integer-labelled')
        elif name.endswith('_surface'):
            mesh = read_legacy_vtk(path, input_frame='RAS')
            meshes[name] = {'vertices': len(mesh.vertices), 'triangles': len(mesh.faces)}
    return {
        'status': 'passed', 'scope': 'Python execution and output integrity; source agreement unverified',
        'manifest': str(manifest_path), 'manifest_sha256': sha256(manifest_path),
        'stage_count': len(expected), 'verified_file_count': len(paths),
        'finite_ras_nifti_count': len(images),
        'nonfinite_intermediate_images_requiring_review': nonfinite_diagnostic_images,
        'meshes': meshes, 'products': {name: str(path) for name, path in products.items()},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    try:
        report = audit(args.run)
        code = 0
    except (OSError, ValueError, KeyError, TypeError) as exc:
        report = {'status': 'failed', 'error': str(exc)}
        code = 2
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return code


if __name__ == '__main__':
    raise SystemExit(main())
