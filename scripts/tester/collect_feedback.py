#!/usr/bin/env python3
"""Collect BASE manifests, logs and version details; exclude image/surface files."""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import zipfile

PACKAGES = ('pymacruise','numpy','scipy','nibabel','numba','llvmlite','antspyx')
ROOT_FILES = ('run.json','checks.json','tester-run-summary.json','comparison.json','tester-run.log',
              'image-diagnostics.json','review-geometry.json')


def image_diagnostics(run, manifest):
    """Return aggregate geometry/value fingerprints, with no voxel arrays or samples."""
    import nibabel as nib
    import numpy as np
    selected={}
    for stage in ('n4','base-preprocessing'):
        for item in manifest.get('stages',{}).get(stage,{}).get('outputs',[]):
            if item['path'].endswith('.nii.gz'):
                selected[stage+'/'+Path(item['path']).name]=Path(item['path'])
    selected.update({'input/'+key:Path(item['path']) for key,item in manifest.get('inputs',{}).items()
                     if key in ('t1','brain_mask','maseg')})
    selected.update({'final/'+Path(item['path']).name:Path(item['path'])
                     for item in manifest.get('final_outputs',[]) if Path(item['path']).name=='SEG_refine.nii.gz'})
    reports={}
    for name,path in selected.items():
        try:
            im=nib.load(path);data=np.asanyarray(im.dataobj)
            if not np.isfinite(data).all():raise ValueError('Image contains nonfinite values')
            canonical=np.ascontiguousarray(data,dtype=data.dtype.newbyteorder('<'))
            flat=np.asarray(data,dtype=np.float64)
            reports[name]={'shape':list(im.shape),'dtype':str(data.dtype),
                'voxel_digest_encoding':'C-contiguous little-endian stored values after NIfTI scaling',
                'voxel_sha256':hashlib.sha256(memoryview(canonical).cast('B')).hexdigest(),
                'minimum':float(flat.min()),'maximum':float(flat.max()),'mean':float(flat.mean()),
                'standard_deviation':float(flat.std()),'nonzero_voxels':int(np.count_nonzero(data)),
                'affine':im.affine.tolist(),'qform':im.get_qform().tolist(),'sform':im.get_sform().tolist(),
                'qform_code':int(im.header['qform_code']),'sform_code':int(im.header['sform_code'])}
        except Exception as error:reports[name]={'error':str(error)}
    return {'scope':'Aggregate image statistics and geometry only; no image arrays or voxel samples',
            'images':reports}


def environment():
    versions = {}
    for name in PACKAGES:
        try: versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: versions[name] = None
    backend=[]
    try:
        distribution=importlib.metadata.distribution('antspyx')
        for item in distribution.files or ():
            if str(item).endswith(('.so','.pyd','.dylib')):
                path=distribution.locate_file(item);digest=hashlib.sha256()
                with path.open('rb') as stream:
                    for block in iter(lambda:stream.read(1048576),b''):digest.update(block)
                backend.append({'distribution_path':str(item),'sha256':digest.hexdigest()})
    except (importlib.metadata.PackageNotFoundError,OSError):pass
    return {'python':platform.python_version(),'executable':sys.executable,
            'platform':platform.platform(),'machine':platform.machine(),'packages':versions,
            'ants_native_backend_files':backend}


def collect(run: Path, *, t1: Path | None = None, verify: bool = True, bundle: Path | None = None):
    run = run.expanduser().resolve()
    if not run.is_dir(): raise ValueError(f'Run directory does not exist: {run}')
    bundle = bundle or Path(__file__).resolve().parent
    verification = {'performed':False,'passed':False}
    if verify:
        try:
            manifest = json.loads((run/'run.json').read_text())
            reference = t1 or Path(manifest['inputs']['t1']['path'])
            result = subprocess.run([sys.executable,str(bundle/'verify-base-run.py'),
                '--run',str(run),'--t1',str(reference),'--report',str(run/'checks.json')],
                capture_output=True,text=True,timeout=180)
            verification = {'performed':True,'passed':result.returncode==0,'returncode':result.returncode}
            if not (run/'checks.json').exists():
                (run/'checks.json').write_text(json.dumps({'artifact_and_grid_checks_passed':False,
                    'error':(result.stderr or result.stdout)[-6000:]},indent=2)+'\n')
        except Exception as error:
            verification['error'] = str(error)
            (run/'checks.json').write_text(json.dumps({'artifact_and_grid_checks_passed':False,
                'error':str(error)},indent=2)+'\n')
    try:
        manifest=json.loads((run/'run.json').read_text())
        if manifest.get('inputs'):
            if (run/'image-diagnostics.json').is_symlink():
                raise ValueError('Refusing to write image diagnostics through a symlink')
            diagnostics=image_diagnostics(run,manifest)
            (run/'image-diagnostics.json').write_text(json.dumps(diagnostics,indent=2,allow_nan=False)+'\n')
    except Exception as error:
        if verify and not (run/'image-diagnostics.json').is_symlink():
            (run/'image-diagnostics.json').write_text(json.dumps({'unavailable':str(error)},indent=2)+'\n')
    summary_path=run/'tester-run-summary.json'
    if verify and summary_path.is_file() and not summary_path.is_symlink():
        try:
            execution=json.loads(summary_path.read_text())
            pipeline_code=execution.get('pipeline_returncode',execution['returncode'])
            execution.update(pipeline_returncode=pipeline_code,
                verification_passed=verification['passed'],
                returncode=pipeline_code or (0 if verification['passed'] else 1))
            summary_path.write_text(json.dumps(execution,indent=2)+'\n')
        except (ValueError,KeyError,TypeError):pass
    summary = {'environment':environment(),'verification':verification,
        'included':'run/check/comparison JSON, aggregate image diagnostics, review geometry, launcher and stage logs, stage provenance',
        'excluded':'NIfTI volumes, meshes, screenshots, raw arrays and input images',
        'metadata_notice':'Local file paths and input/output hashes remain in manifests and logs.'}
    archive = run/(run.name+'-feedback.zip')
    paths = [run/name for name in ROOT_FILES]
    paths += list((run/'stages').glob('*/attempt-*/stage.log'))
    paths += list((run/'stages').glob('*/attempt-*/provenance.json'))
    included = []
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as stream:
        for path in sorted(set(paths)):
            if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(run):continue
            name = path.relative_to(run).as_posix()
            stream.write(path,name);included.append(name)
        summary['files'] = included
        stream.writestr('feedback-environment.json',json.dumps(summary,indent=2)+'\n')
    return archive, summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',required=True,type=Path)
    parser.add_argument('--t1',type=Path,help='Original T1; otherwise read its path from run.json')
    args = parser.parse_args(argv)
    archive,summary = collect(args.run,t1=args.t1)
    print('Feedback archive:',archive)
    print('Run checks:', 'passed' if summary['verification']['passed'] else 'failed or incomplete; diagnostics included')
    return 0

if __name__=='__main__':raise SystemExit(main())
