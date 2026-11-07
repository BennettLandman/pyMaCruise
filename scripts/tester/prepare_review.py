#!/usr/bin/env python3
"""Create paired review copies on the original T1 grid and in physical RAS mm."""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import nibabel as nib
import numpy as np
from pymacruise.surfaces import SurfaceMesh, read_legacy_vtk, write_legacy_vtk


def geometry_module():
    folder=Path(__file__).resolve().parent
    path=folder/'compare_base.py'
    if not path.is_file():path=folder/'compare-with-macruise.py'
    spec=importlib.util.spec_from_file_location('review_comparison',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def world_vertices(vertices, affine):
    result=np.einsum('ij,nj->ni',affine[:3,:3],vertices,optimize=False)+affine[:3,3]
    if not np.isfinite(result).all():raise ValueError('Surface transform produced nonfinite coordinates')
    return result


def prepare(run, source, output, *, surfaces=True):
    run=run.expanduser().resolve();source=source.expanduser().resolve();output=output.expanduser().resolve()
    manifest=json.loads((run/'run.json').read_text())
    if manifest['status']!='complete':raise ValueError('Complete the run before preparing paired review outputs')
    if output.exists():raise FileExistsError('Use a new review output directory')
    if (run/'review-geometry.json').is_symlink():raise ValueError('Refusing to write review geometry through a symlink')
    products={Path(item['path']).name:Path(item['path']) for item in manifest['final_outputs']}
    t1=nib.load(manifest['inputs']['t1']['path']);raw=nib.load(source/'Raw/Seg.nii.gz')
    if not np.array_equal(np.asarray(nib.load(source/'Raw/T1.nii.gz').dataobj,dtype=np.float32),np.asarray(t1.dataobj,dtype=np.float32)):
        raise ValueError('Original Raw/T1 voxel values differ from this run input')
    if not np.array_equal(np.asarray(raw.dataobj),np.asarray(nib.load(manifest['inputs']['maseg']['path']).dataobj)):
        raise ValueError('Original Raw/Seg voxel values differ from this run input MaSeg')
    original=nib.load(source/'SegRefine/SEG_refine.nii.gz');candidate=nib.load(products['SEG_refine.nii.gz'])
    geometry=geometry_module()
    report={'purpose':'Paired geometry review; no reconstruction, registration, smoothing or label edits',
        'source_run_converged':manifest.get('converged'),
        'original_vs_python_geometry':geometry.geometry_metrics(original,candidate),
        'label_sampling':'none; voxel arrays preserved; only sub-0.0001 mm header rounding is standardized',
        'surface_frame':'physical RAS millimetres','fitted_registration':False,'surfaces':{}}
    for name,image in [('original',original),('python',candidate)]:
        check=geometry.geometry_metrics(image,t1)
        if not check['same_lattice_within_0_0001_mm']:
            raise ValueError(name+' SegRefine differs from the input grid by more than rounding; refusing header-only copying')
        values=np.asanyarray(image.dataobj)
        if not np.isfinite(values).all() or np.any(values<0) or not np.equal(values,np.rint(values)).all():
            raise ValueError(name+' SegRefine contains invalid labels')
        if values.max()>np.iinfo(np.int32).max:raise ValueError('Labels exceed int32 storage')
    # Validate the complete source coordinate convention before writing copies.
    meshes={}
    if surfaces:
        for name in ('inner','central','outer'):
            old=read_legacy_vtk(source/f'Surfaces/target_image_GMimg_{name}Surf.vtk',input_frame='RAS')
            world=world_vertices(old.vertices,raw.affine)
            with (source/f'Surfaces_FreeView/target_image_GMimg_{name}Surf.asc').open() as stream:
                next(stream);count,_=map(int,next(stream).split())
                asc=np.array([list(map(float,next(stream).split()[:3])) for _ in range(count)])
            if asc.shape!=world.shape:raise ValueError('Original raw VTK and ASC vertex counts differ')
            error=float(np.linalg.norm(asc-world,axis=1).max())
            if not np.isfinite(error) or error>=2e-5:
                raise ValueError('Original raw VTK does not map to its ASC physical export')
            meshes['original_'+name]=SurfaceMesh(world,old.faces,point_data=old.point_data)
            meshes['python_'+name]=read_legacy_vtk(products[name+'_surface.vtk'])
            report['surfaces'][name]={'original_voxel_to_ASC_max_error_mm':error,
                'faces_preserved':True,'original_vertices':len(world),
                'python_vertices':len(meshes['python_'+name].vertices)}
    output.mkdir(parents=True,exist_ok=False)
    for name,image in [('original',original),('python',candidate)]:
        values=np.asarray(image.dataobj,dtype=np.int32)
        header=t1.header.copy();header.set_data_dtype(np.int32);header.set_slope_inter(None,None)
        kind=nib.Nifti2Image if isinstance(t1,nib.Nifti2Image) else nib.Nifti1Image
        result=kind(values,t1.affine,header=header)
        result.set_qform(t1.get_qform(),code=int(t1.header['qform_code']))
        result.set_sform(t1.get_sform(),code=int(t1.header['sform_code']))
        path=output/(name+'_SEG_refine.nii.gz');nib.save(result,path)
        saved=nib.load(path)
        np.testing.assert_array_equal(np.asarray(saved.dataobj),values)
        np.testing.assert_array_equal(saved.affine,t1.affine)
    for name,mesh in meshes.items():write_legacy_vtk(mesh,output/(name+'_RAS_mm.vtk'))
    report['review_output_directory']=str(output)
    text=json.dumps(report,indent=2,allow_nan=False)+'\n'
    (output/'review-geometry.json').write_text(text)
    (run/'review-geometry.json').write_text(text)
    (output/'README.txt').write_text(
        'Load original_SEG_refine.nii.gz and python_SEG_refine.nii.gz over the original input T1.\n'
        'Use nearest-neighbor label display. Both copies preserve every original label value.\n'
        'Load original_<inner|central|outer>_RAS_mm.vtk and python_<inner|central|outer>_RAS_mm.vtk together.\n'
        'All paired meshes use physical RAS millimetres. Original raw VTK voxel indices were converted\n'
        'through their reference affine and checked against the original ASC export.\n'
        'These copies standardize geometry for review; anatomical differences and nonconvergence remain.\n')
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',required=True,type=Path)
    parser.add_argument('--original-output',required=True,type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--labels-only',action='store_true')
    args=parser.parse_args(argv)
    try:report=prepare(args.run,args.original_output,args.output_dir,surfaces=not args.labels_only)
    except Exception as error:parser.exit(1,f'Review preparation failed: {error}\n')
    print('Paired review copies:',args.output_dir)
    print('Source run converged:',report['source_run_converged'])
    print('Original vs Python maximum corner displacement (mm):',report['original_vs_python_geometry']['maximum_corner_displacement_mm'])
    return 0


if __name__=='__main__':raise SystemExit(main())
