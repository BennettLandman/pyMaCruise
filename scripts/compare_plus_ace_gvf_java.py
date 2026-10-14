#!/usr/bin/env python3
"""Compare held-fixed saved plus ACE and GVF stages with selected Java source."""
from __future__ import annotations
import argparse,json,re,struct
from pathlib import Path
import numpy as np

def read_volumes(path,shape,count):
    data=path.read_bytes();offset=0;out=[]
    for _ in range(count):
        found=struct.unpack('>III',data[offset:offset+12]);offset+=12
        if found!=shape:raise SystemExit(f'{path.name}: output shape {found} != {shape}')
        n=int(np.prod(shape));a=np.frombuffer(data,dtype='>f4',count=n,offset=offset).astype(np.float32);offset+=n*4;out.append(a.reshape(shape))
    if offset!=len(data):raise SystemExit(f'{path.name}: trailing bytes')
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reference-dir',type=Path,required=True);p.add_argument('--report',type=Path,required=True);a=p.parse_args()
    root=a.reference_dir;provenance=json.loads((root/'input-provenance.json').read_text());shape=tuple(provenance['shape_xyz'])
    java_ace=read_volumes(root/'ace-output.bin',shape,3)
    ace_names=('ace-enhanced','ace-skeleton','ace-thinned-skeleton')
    ace_rows={}
    for name,java in zip(ace_names,java_ace,strict=True):
        py=np.load(root/f'python-{name}.npy').astype(np.float32);err=np.abs(java.astype(np.float64)-py.astype(np.float64))
        ace_rows[name]={'exact':bool(np.array_equal(java,py)),'changed_voxels':int(np.count_nonzero(java!=py)),'voxel_count':int(java.size),'maximum_absolute_error':float(err.max(initial=0)),'mean_absolute_error':float(err.mean())}
    raw=(root/'gvf-output.bin').read_bytes();found=struct.unpack('>III',raw[:12])
    if found!=shape:raise SystemExit(f'GVF shape mismatch: {found} != {shape}')
    java_gvf=np.frombuffer(raw,dtype='>f4',offset=12).astype(np.float32)
    if java_gvf.size!=int(np.prod(shape))*3:raise SystemExit('GVF output has an unexpected payload length')
    java_gvf=java_gvf.reshape((*shape,3));python_gvf=np.load(root/'python-gvf-vector.npy').astype(np.float32)
    err=np.abs(java_gvf.astype(np.float64)-python_gvf.astype(np.float64))
    gvf={'exact':bool(np.array_equal(java_gvf,python_gvf)),'changed_values':int(np.count_nonzero(java_gvf!=python_gvf)),'value_count':int(java_gvf.size),'maximum_absolute_error':float(err.max(initial=0)),'mean_absolute_error':float(err.mean()),'accepted_by_existing_reference_limits':bool(float(err.mean())<=1e-5 and float(err.max(initial=0))<=5e-5)}
    def peak(name):
        text=(root/name).read_text();matches=re.findall(r"peak_sampled_rss_gb=([0-9.]+)",text)
        if not matches:raise SystemExit(f'{name}: watchdog peak RSS missing')
        value=float(matches[-1]);return value if value>=0.01 else None
    report={'synthetic_only':provenance['synthetic_only'],'evidence_category':provenance['evidence_category'],'comparison':'selected MaCRUISEplus Java ACE mode 3 and VandyCRUISE GVF vs current plus workflow on exact held-fixed inputs','scope':'saved Python TGDM outputs are hash-checked; ACE and GVF inputs are reconstructed from the same crop inputs and plus settings','workflow_run_signature':provenance['workflow_run_signature'],'implementation_signature':provenance['implementation_signature'],'shape_xyz':list(shape),'coordinate_system':'RAS','parameters':provenance['ace_gvf_parameters'],'java_peak_sampled_rss_gb':{'ace':peak('ace-java-watchdog.log'),'gvf':peak('gvf-java-watchdog.log')},'java_rss_limit_gb':8,'ace_outputs':ace_rows,'gvf_output':gvf,'passed':all(v['exact'] for v in ace_rows.values()) and gvf['accepted_by_existing_reference_limits']}
    a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,indent=2))
    if not report['passed']:raise SystemExit('plus ACE/GVF outputs missed selected-source comparison tolerances')
if __name__=='__main__':main()
