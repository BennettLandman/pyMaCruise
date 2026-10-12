#!/usr/bin/env python3
"""Compare Java CropBasedOnGM on legacy-sized transfer maps to Python physical transfer + crop."""
from __future__ import annotations
import argparse, hashlib, json, os, subprocess
from pathlib import Path
import nibabel as nib
import numpy as np
from compare_human_jist_base_crop import (
    RUNTIME_JAR, WRAPPER_SOURCE, classpath, jar_entry_sha256,
)
from pymacruise.cropping import find_crop_parameters, crop_image
from pymacruise.spatial_transfer import transfer_to_reference_space

ROOT=Path(__file__).resolve().parents[1]
HARNESS=ROOT/'tests/java_reference/P41zUnequalGridTransferCropReferenceHarness.java'
TARGETS=('csf','gm','wm_fill','wm_mask')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write_raw(path,values,integer=False):
 dtype='>i4' if integer else '>f4'
 np.asarray(values,dtype=dtype).ravel(order='C').tofile(path)
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--reference',type=Path,required=True)
 p.add_argument('--atlas-gm',type=Path,required=True)
 p.add_argument('--dura-dir',type=Path,required=True,help='directory containing lesion_toads_{csf,gm,wm_fill,wm_mask}.nii.gz')
 p.add_argument('--edge-width',type=int,default=5)
 p.add_argument('--output-dir',type=Path,required=True)
 a=p.parse_args(); refp=a.reference.resolve(); gmp=a.atlas_gm.resolve(); duradir=a.dura_dir.resolve(); out=a.output_dir.resolve()
 ref=nib.load(str(refp)); gm=nib.load(str(gmp))
 if ref.shape!=gm.shape or not np.allclose(ref.affine,gm.affine,rtol=0,atol=1e-5): raise ValueError('T1 reference and Atlas_GM must share the source reference grid')
 if nib.aff2axcodes(ref.affine)!=('R','A','S'): raise ValueError('reference must be RAS')
 sources={n:duradir/f'lesion_toads_{n}.nii.gz' for n in TARGETS}
 imgs={n:nib.load(str(path)) for n,path in sources.items()}
 shapes={tuple(i.shape) for i in imgs.values()}
 if len(shapes)!=1: raise ValueError('Dura source maps must share a grid')
 source_shape=next(iter(shapes))
 if any(nib.aff2axcodes(i.affine)!=('R','A','S') for i in imgs.values()): raise ValueError('source maps must be RAS')
 out.mkdir(parents=True,exist_ok=True); inp=out/'java-input'; jout=out/'java-output'; classes=out/'classes'; home=out/'java-home'
 for d in (inp,jout,classes,home):d.mkdir(parents=True,exist_ok=True)
 gm_data=np.asarray(gm.dataobj)
 write_raw(inp/'atlas_gm.i32be',gm_data.astype(np.int32),True)
 for name,img in imgs.items(): write_raw(inp/f'{name}.f32be',np.asarray(img.dataobj,dtype=np.float32))
 cp=classpath()
 subprocess.run(['javac','-classpath',cp,'-d',str(classes),str(HARNESS)],check=True)
 java_cp=str(classes)+os.pathsep+cp
 with (out/'java-crop.log').open('w',encoding='utf-8') as log:
  subprocess.run(['java','-Xmx2g','-Djava.awt.headless=true',f'-Duser.home={home}','-cp',java_cp,'P41zUnequalGridTransferCropReferenceHarness',str(inp),str(jout),*(str(v) for v in gm.shape),*(str(v) for v in source_shape),str(a.edge_width)],check=True,stdout=log,stderr=subprocess.STDOUT)
 crop=find_crop_parameters(gm,threshold=0.0,edge_width=a.edge_width)
 rows=[]
 for name,img in imgs.items():
  java=np.fromfile(jout/f'{name}.f32be',dtype='>f4').astype(np.float32)
  jshape=tuple(int(x) for x in (jout/'bounds.txt').read_text().strip().split(';')[-1].split(','))
  java=java.reshape(jshape,order='C')
  python_full=transfer_to_reference_space(img,ref,cval=0.0)
  python_crop=crop_image(python_full,crop)
  pv=np.asarray(python_crop.dataobj,dtype=np.float32)
  delta=np.abs(java-pv)
  rows.append({'map':name,'input_shape':list(img.shape),'crop_shape':list(jshape),'java_vs_python_physical_transfer_then_crop_differing_voxels':int(np.count_nonzero(java!=pv)),'voxel_count':int(java.size),'max_abs_difference':float(delta.max(initial=0.0)),'mean_abs_difference':float(delta.mean(dtype=np.float64)),'java_nonzero_voxels':int(np.count_nonzero(java)),'python_nonzero_voxels':int(np.count_nonzero(pv)),'source_sha256':sha(sources[name])})
 report={'case':'P41z iteration-11 synthetic transfer followed by selected Java crop','evidence_category':'synthetic selected Java crop comparison; no full workflow equivalence claim','java_cropper_class_sha256':jar_entry_sha256(RUNTIME_JAR,'edu/jhmi/rad/medic/utilities/CubicVolumeCropper.class'),'crop_parameters_class_sha256':jar_entry_sha256(RUNTIME_JAR,'edu/jhmi/rad/medic/utilities/CropParameters.class'),'selected_wrapper':str(WRAPPER_SOURCE),'selected_wrapper_sha256':sha(WRAPPER_SOURCE),'reference_shape_xyz':list(ref.shape),'dura_source_shape_xyz':list(source_shape),'python_crop_bounds_inclusive':[list(x) for x in crop.bounds],'java_crop_bounds_and_shape':(jout/'bounds.txt').read_text().strip(),'atlas_gm_sha256':sha(gmp),'reference_sha256':sha(refp),'dura_source_geometry':list(nib.aff2axcodes(imgs['gm'].affine)),'results':rows,'all_maps_exact':all(x['java_vs_python_physical_transfer_then_crop_differing_voxels']==0 for x in rows)}
 report_path=out/'comparison.json';report_path.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,indent=2));print('report='+str(report_path))
 return 0
if __name__=='__main__': raise SystemExit(main())
