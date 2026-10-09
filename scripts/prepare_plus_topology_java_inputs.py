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
from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_tissue_update_stage, run_plus_topology_stage
from pymacruise.topology import _erode_multi_atlas_wm


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workflow-root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args(); root=args.workflow_root.resolve()
    manifest=json.loads((root/'run.json').read_text())
    records=manifest['stages']['transfer_crop']['outputs']
    def load(name): return nib.load(records[name]['path'])
    settings=dict(manifest['settings']['plus_pipeline'])
    settings['boundary_overlook_labels']=tuple(settings['boundary_overlook_labels'])
    options=PlusPipelineOptions(**settings)
    tissue=run_plus_tissue_update_stage(load('multi_atlas_gm'),load('multi_atlas_wm'),load('multi_atlas_mask'),load('toads_gm'),load('toads_wm_fill'),load('toads_wm_mask'))
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
    report={
        'synthetic_only':True,'workflow_manifest':str(root/'run.json'),'workflow_manifest_sha256':sha256(root/'run.json'),
        'workflow_run_signature':manifest['run_signature'],'shape_xyz':list(shape),'coordinate_system':'RAS',
        'affine_ras':tissue.wm.affine.tolist(),'settings':options.as_dict(),
        'held_fixed_input_binary_sha256':sha256(input_path),
        'crop_inputs':{name:{'path':records[name]['path'],'sha256':sha256(Path(records[name]['path']))} for name in ('multi_atlas_gm','multi_atlas_wm','multi_atlas_mask','toads_gm','toads_wm_fill','toads_wm_mask')},
        'outputs':['corrected','mawm_eroded','atlas_refined'],
        'scope':'selected Java topology correction compares the corrected scalar image; the MaWM erosion and Refine_Topo_by_MA post-step are separately compared using the preserved MASImorph source.'}
    (args.output_dir/'topology-provenance.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
