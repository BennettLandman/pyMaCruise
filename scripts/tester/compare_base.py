#!/usr/bin/env python3
"""Compare a completed BASE run with the same case's selected MaCRUISE outputs."""
from __future__ import annotations
import argparse
import hashlib
import json
import itertools
from pathlib import Path
import nibabel as nib
import numpy as np


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):digest.update(block)
    return digest.hexdigest()


def label_metrics(original,candidate):
    a=np.asarray(original);b=np.asarray(candidate)
    if a.shape!=b.shape:raise ValueError('Refined segmentations have different voxel dimensions')
    for values in (a,b):
        if not np.isfinite(values).all() or np.any(values<0) or not np.equal(values,np.rint(values)).all():
            raise ValueError('Segmentation labels must be finite, nonnegative integers')
    labels=np.union1d(np.unique(a),np.unique(b))
    ia=np.searchsorted(labels,a.ravel());ib=np.searchsorted(labels,b.ravel());n=len(labels)
    ca=np.bincount(ia,minlength=n);cb=np.bincount(ib,minlength=n)
    same=ia==ib;intersection=np.bincount(ia[same],minlength=n)
    per={str(int(k)):{'dice':float(2*intersection[j]/(ca[j]+cb[j])),
        'original_voxels':int(ca[j]),'python_voxels':int(cb[j])}
        for j,k in enumerate(labels) if k!=0}
    bad=[k for k,v in per.items() if v['dice']<.97]
    return {'differing_voxels':int(np.count_nonzero(a!=b)),'voxel_count':int(a.size),
        'mean_label_dice':float(np.mean([v['dice'] for v in per.values()])) if per else 1.,
        'labels_below_0_97':bad,'per_label':per,'dice_gate_passed':not bad}


def geometry_metrics(original, candidate):
    """Quantify physical displacement rather than treating rounding as misalignment."""
    corners=np.array(list(itertools.product(*[(0,n-1) for n in candidate.shape])),dtype=float)
    delta=original.affine-candidate.affine
    shifts=np.einsum('ij,nj->ni',delta[:3,:3],corners,optimize=False)+delta[:3,3]
    maximum=float(np.linalg.norm(shifts,axis=1).max())
    return {'original_shape':list(original.shape),'python_shape':list(candidate.shape),
        'original_affine':original.affine.tolist(),'python_affine':candidate.affine.tolist(),
        'original_qform_code':int(original.header['qform_code']),
        'python_qform_code':int(candidate.header['qform_code']),
        'original_sform_code':int(original.header['sform_code']),
        'python_sform_code':int(candidate.header['sform_code']),
        'affines_bit_identical':bool(np.array_equal(original.affine,candidate.affine)),
        'maximum_corner_displacement_mm':maximum,
        'same_lattice_within_0_0001_mm':original.shape==candidate.shape and maximum<=1e-4}


def compare(run,source,*,surfaces=False,report_path=None):
    run=run.resolve();source=source.resolve();report_path=report_path or run/'comparison.json'
    manifest=json.loads((run/'run.json').read_text())
    if manifest['status']!='complete':raise ValueError('Complete the Python run before comparing anatomy')
    files={Path(r['path']).name:Path(r['path']) for r in manifest['final_outputs']}
    reference=nib.load(manifest['inputs']['t1']['path']);raw_seg=nib.load(source/'Raw/Seg.nii.gz')
    if not np.array_equal(np.asarray(raw_seg.dataobj),np.asarray(nib.load(manifest['inputs']['maseg']['path']).dataobj)):
        raise ValueError('Original Raw/Seg does not match this run’s input MaSeg voxel array')
    if not np.array_equal(np.asarray(nib.load(source/'Raw/T1.nii.gz').dataobj,dtype=np.float32),np.asarray(reference.dataobj,dtype=np.float32)):
        raise ValueError('Original Raw/T1 does not match this run’s input T1 voxel array')
    py=nib.load(files['SEG_refine.nii.gz'])
    if py.shape!=reference.shape or not np.array_equal(py.affine,reference.affine):
        raise ValueError('Python final refined segmentation is not on the original T1 grid')
    original=nib.load(source/'SegRefine/SEG_refine.nii.gz')
    data={'comparison_status':'running','scope':'same-case original-versus-Python comparison; no fitted registration or cleanup',
        'python_run':str(run),'original_output':str(source),'manifest_sha256':sha(run/'run.json'),
        't1_and_maseg_voxel_identity_verified':True,'segmentation':label_metrics(original.dataobj,py.dataobj),
        'segmentation_geometry':geometry_metrics(original,py),
        'surfaces':{},'project_gates':{'label_dice_minimum':.97,'surface_symmetric_mean_maximum_mm':.01,'sampled_maximum_mm':1.},
        'limitations':['Software checks do not establish scientific or clinical equivalence.',
            'Surface distances sample every vertex against the opposite triangle mesh; the maximum is not a continuous Hausdorff bound.']}
    mask_paths=list((source/'ToadsSeperate/MultiAtlasStripping').glob('*_MAstrip_mask.nii.gz'))
    mask_identity={'checked':False,'matched':None,'reason':'Original mask not identified automatically; confirm matching skull-stripping inputs.'}
    if len(mask_paths)==1:
        original_mask=np.asarray(nib.load(mask_paths[0]).dataobj)>0
        input_mask=np.asarray(nib.load(manifest['inputs']['brain_mask']['path']).dataobj)>0
        matched=original_mask.shape==input_mask.shape and np.array_equal(original_mask,input_mask)
        mask_identity={'checked':True,'matched':matched,'source':'original MultiAtlas stripping mask',
            'different_voxels':int(np.count_nonzero(original_mask!=input_mask)) if original_mask.shape==input_mask.shape else None}
    data['original_brain_mask_identity']=mask_identity
    if mask_identity['matched'] is not True:
        data['limitations'].append('Original brain-mask identity is unconfirmed or differs; masking may confound the comparison.')
    def save():
        temp=report_path.with_suffix('.tmp');temp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');temp.replace(report_path)
    save()
    if surfaces:
        from pymacruise.jist_intersector import JistSurfaceTree
        from pymacruise.surfaces import read_legacy_vtk
        from pymacruise.validation import mesh_topology_metrics
        for name in ('inner','central','outer'):
            print('Comparing',name,'surface: all vertices against opposite triangles.',flush=True)
            old=read_legacy_vtk(source/f'Surfaces/target_image_GMimg_{name}Surf.vtk',input_frame='RAS')
            world=np.einsum('ij,nj->ni',raw_seg.affine[:3,:3],old.vertices,optimize=False)+raw_seg.affine[:3,3]
            with (source/f'Surfaces_FreeView/target_image_GMimg_{name}Surf.asc').open() as stream:
                next(stream);count,_=map(int,next(stream).split())
                asc=np.array([list(map(float,next(stream).split()[:3])) for _ in range(count)])
            if asc.shape!=world.shape:raise ValueError('Original raw VTK and FreeView ASC have different vertex counts')
            frame_error=float(np.linalg.norm(asc-world,axis=1).max())
            if not np.isfinite(frame_error) or frame_error>=2e-5:raise ValueError('Original VTK voxel-to-RAS mapping does not match the supplied ASC; inspect this source version’s coordinate convention')
            new=read_legacy_vtk(files[name+'_surface.vtk'])
            tree=JistSurfaceTree(new.vertices,new.faces);_,da,_=tree.query_many(world);del tree
            tree=JistSurfaceTree(world,old.faces);_,db,_=tree.query_many(new.vertices);del tree
            mean=float((da.mean()+db.mean())/2);maximum=float(max(da.max(),db.max()))
            data['surfaces'][name]={'symmetric_mean_vertex_to_triangle_mm':mean,'sampled_hausdorff_mm':maximum,
                'original_frame_to_asc_max_error_mm':frame_error,'distance_gate_passed':mean<=.01 and maximum<=1,
                'original_topology':mesh_topology_metrics(world,old.faces),'python_topology':mesh_topology_metrics(new.vertices,new.faces)}
            save()
    data['comparison_status']='complete'
    data['segmentation_and_surface_gates_passed']=(data['segmentation']['dice_gate_passed'] and all(v['distance_gate_passed'] for v in data['surfaces'].values())) if surfaces else None
    save();return data


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',required=True,type=Path)
    parser.add_argument('--original-output',required=True,type=Path,help='Original MaCRUISE Output folder containing Raw, SegRefine and Surfaces')
    parser.add_argument('--surfaces',action='store_true',help='Also compute all-vertex surface distances; takes additional time')
    args=parser.parse_args(argv)
    try:
        data=compare(args.run,args.original_output,surfaces=args.surfaces)
    except Exception as error:
        if args.run.is_dir():(args.run/'comparison.json').write_text(json.dumps({'comparison_status':'failed','error':str(error)},indent=2)+'\n')
        parser.exit(1,f'Comparison failed: {error}\n')
    print('Comparison saved:',args.run/'comparison.json')
    print('Differing refined-label voxels:',data['segmentation']['differing_voxels'])
    print('Label Dice gate:',data['segmentation']['dice_gate_passed'])
    if args.surfaces:print('Combined project gates:',data['segmentation_and_surface_gates_passed'])
    print('A completed comparison can report failed agreement gates; inspect the measurements.')
    return 0

if __name__=='__main__':raise SystemExit(main())
