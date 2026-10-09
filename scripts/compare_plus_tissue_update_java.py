#!/usr/bin/env python3
"""Compare preserved CombineMAandTOADS Java source with the plus Python stage."""
from __future__ import annotations
import argparse,gzip,json
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reference-dir',type=Path,required=True);p.add_argument('--report',type=Path,required=True);a=p.parse_args()
    root=a.reference_dir;provenance=json.loads((root/'topology-provenance.json').read_text());shape=tuple(provenance['shape_xyz'])
    rows={}
    for name in ('gm','wm','mask'):
        with gzip.open(root/f'java-{name}.raw.gz','rb') as f:
            header=tuple(int.from_bytes(f.read(4),'big') for _ in range(3));java=np.frombuffer(f.read(),dtype='>f4').astype(np.float32)
        if header!=shape or java.size!=int(np.prod(shape)):raise SystemExit(f'{name}: Java output shape mismatch')
        java=java.reshape(shape,order='C');python=np.load(root/f'python-tissue-{name}.npy').astype(np.float32)
        diff=np.abs(java.astype(np.float64)-python.astype(np.float64))
        rows[name]={'exact':bool(np.array_equal(java,python)),'changed_voxels':int(np.count_nonzero(java!=python)),'voxel_count':int(java.size),'maximum_absolute_error':float(diff.max(initial=0)),'mean_absolute_error':float(diff.mean())}
    report={'synthetic_only':True,'comparison':'preserved selected Java CombineMAandTOADS.execute() vs Python combine_multi_atlas_toads on workflow crop inputs','scope':'held-fixed six atlas/TOADS input images from the saved plus run; selected default ErodeMethod=3D and GMconfident=0.001','workflow_run_signature':provenance['workflow_run_signature'],'shape_xyz':list(shape),'coordinate_system':'RAS','outputs':rows,'all_outputs_exact':all(r['exact'] for r in rows.values())}
    a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,indent=2))
    if not report['all_outputs_exact']:raise SystemExit('Java and Python plus tissue products differ')
if __name__=='__main__':main()
