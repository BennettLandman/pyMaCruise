#!/usr/bin/env python3
"""Prepare exact crop-grid inputs for the selected plus Java topology comparison."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import struct
import nibabel as nib
import numpy as np
from pymacruise.cropping import CropParameters, crop_image
from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_tissue_update_stage, run_plus_topology_stage
from pymacruise.topology import _erode_multi_atlas_wm


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workflow-root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument(
        '--multi-atlas-brain-mask', type=Path,
        help='Optional full-grid Atlas_BrainTissueMask; crop with saved bounds instead of using WM-only Atlas_Mask.',
    )
    args=parser.parse_args(); root=args.workflow_root.resolve()
    manifest=json.loads((root/'run.json').read_text())
    records=manifest['stages']['transfer_crop']['outputs']
    def load(name): return nib.load(records[name]['path'])
    mask_image = load('multi_atlas_mask')
    mask_source = records['multi_atlas_mask']['path']
    if args.multi_atlas_brain_mask is not None:
        crop = manifest['stages']['transfer_crop']['crop_parameters']
        parameters = CropParameters(
            original_shape=tuple(crop['original_shape']),
            bounds=tuple(tuple(pair) for pair in crop['bounds']),
            background=float(crop['background']),
            original_affine=np.asarray(crop['original_affine'], dtype=np.float64),
        )
        mask_image = crop_image(nib.load(str(args.multi_atlas_brain_mask)), parameters)
        args.output_dir.mkdir(parents=True,exist_ok=True)
        mask_path=args.output_dir/'multi_atlas_brain_mask_cropped.nii.gz'
        nib.save(mask_image,mask_path)
        mask_source=str(mask_path.resolve())
    settings=dict(manifest['settings']['plus_pipeline'])
    settings['boundary_overlook_labels']=tuple(settings['boundary_overlook_labels'])
    options=PlusPipelineOptions(**settings)
    tissue=run_plus_tissue_update_stage(
        load('multi_atlas_gm'), load('multi_atlas_wm'), mask_image,
        load('toads_gm'), load('toads_wm_fill'), load('toads_wm_mask'),
    )
    topology=run_plus_topology_stage(tissue.wm,tissue.mask,load('multi_atlas_wm'),options)
    if topology.atlas_refined is None: raise SystemExit('saved plus settings require MaSeg topology refinement')
    shape=tissue.wm.shape
    arrays=(np.asarray(tissue.wm.dataobj,dtype=np.float32),np.asarray(tissue.mask.dataobj,dtype=np.float32),np.asarray(load('multi_atlas_wm').dataobj,dtype=np.float32))
    args.output_dir.mkdir(parents=True,exist_ok=True)
    input_path=args.output_dir/'topology-arrays.bin'
    with input_path.open('wb') as f:
        f.write(struct.pack('>III',*shape))
        for array in arrays: f.write(np.asarray(array,dtype='>f4').tobytes(order='F'))
    for name,image in (('corrected',topology.corrected),('atlas_refined',topology.atlas_refined)):
        np.save(args.output_dir/f'python-{name}.npy',np.asarray(image.dataobj,dtype=np.float32))
    np.save(args.output_dir/'python-mawm-eroded.npy',_erode_multi_atlas_wm(np.asarray(load('multi_atlas_wm').dataobj)))
    tissue_inputs=(load('multi_atlas_gm'),load('multi_atlas_wm'),mask_image,load('toads_gm'),load('toads_wm_fill'),load('toads_wm_mask'))
    tissue_arrays=tuple(np.asarray(image.dataobj,dtype=np.float32) for image in tissue_inputs)
    tissue_binary=args.output_dir/'tissue-arrays.bin'
    with tissue_binary.open('wb') as f:
        f.write(struct.pack('>III',*shape))
        for array in tissue_arrays: f.write(np.asarray(array,dtype='>f4').tobytes(order='F'))
    for name,image in zip(('gm','wm','mask'),(tissue.gm,tissue.wm,tissue.mask),strict=True):
        np.save(args.output_dir/f'python-tissue-{name}.npy',np.asarray(image.dataobj))
    for name,array in zip(
        ('atlas-gm','atlas-wm','brain-mask','toads-gm','toads-wm-fill','toads-wm-mask'),
        tissue_arrays,
        strict=True,
    ):
        np.save(args.output_dir/f'python-input-{name}.npy',array)
    np.save(args.output_dir/'python-input-maseg.npy',np.asarray(load('segmentation').dataobj))
    report={
        'synthetic_only':True,'workflow_manifest':str(root/'run.json'),'workflow_manifest_sha256':sha256(root/'run.json'),
        'workflow_run_signature':manifest['run_signature'],'shape_xyz':list(shape),'coordinate_system':'RAS',
        'affine_ras':tissue.wm.affine.tolist(),'settings':options.as_dict(),
        'held_fixed_input_binary_sha256':sha256(input_path),'tissue_input_binary_sha256':sha256(tissue_binary),
        'crop_inputs':{
            name: {
                'path': (mask_source if name == 'multi_atlas_mask' else records[name]['path']),
                'sha256': sha256(Path(mask_source if name == 'multi_atlas_mask' else records[name]['path'])),
            }
            for name in ('multi_atlas_gm','multi_atlas_wm','multi_atlas_mask','toads_gm','toads_wm_fill','toads_wm_mask')
        },
        'brain_mask_override_source':(str(args.multi_atlas_brain_mask.resolve()) if args.multi_atlas_brain_mask is not None else None),
        'outputs':['tissue_gm','tissue_wm','tissue_mask','corrected','mawm_eroded','atlas_refined'],
        'scope':'selected Java topology correction compares the corrected scalar image; the MaWM erosion and Refine_Topo_by_MA post-step are separately compared using the preserved MASImorph source.'}
    (args.output_dir/'topology-provenance.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
