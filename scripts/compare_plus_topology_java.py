#!/usr/bin/env python3
"""Compare selected Java plus scalar topology correction with Python outputs."""
from __future__ import annotations
import argparse,gzip,json,re
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reference-dir',type=Path,required=True);p.add_argument('--report',type=Path,required=True);a=p.parse_args()
    root=a.reference_dir; provenance=json.loads((root/'topology-provenance.json').read_text());shape=tuple(provenance['shape_xyz'])
    log=(root/'topology-java-watchdog.log').read_text();m=re.search(r'peak_sampled_rss_gb=([0-9.]+)',log)
    if m is None: raise SystemExit('topology Java watchdog log lacks peak RSS')
    with gzip.open(root/'topology-java.raw.gz','rb') as f:
        header=tuple(int.from_bytes(f.read(4),'big') for _ in range(3))
        java=np.frombuffer(f.read(),dtype='>f4').astype(np.float32)
    if header!=shape or java.size!=int(np.prod(shape)): raise SystemExit('Java output dimensions/payload mismatch')
    java=java.reshape(shape,order='F');python=np.load(root/'python-corrected.npy').astype(np.float32)
    diff=np.abs(java.astype(np.float64)-python.astype(np.float64))
    sampled_rss=float(m.group(1))
    report={
        'synthetic_only':True,'comparison':'selected Java AlgorithmTopologyCorrection vs plus Python topology_correct core','scope':provenance['scope'],
        'workflow_run_signature':provenance['workflow_run_signature'],'shape_xyz':list(shape),'coordinate_system':'RAS','java_peak_sampled_rss_gb':sampled_rss if sampled_rss>=0.01 else None,'java_rss_limit_gb':8,
        'rss_measurement_note':None if sampled_rss>=0.01 else 'The watchdog returned less than 0.01 GB for this short JVM run; treated as below measurement resolution, not a reliable peak.',
        'corrected_scalar_output':{'exact':bool(np.array_equal(java,python)),'changed_voxels':int(np.count_nonzero(java!=python)),'voxel_count':int(java.size),'maximum_absolute_error':float(diff.max(initial=0)),'mean_absolute_error':float(diff.mean()),'label_dice_at_zero':float(2*np.count_nonzero((java<=0)&(python<=0))/(np.count_nonzero(java<=0)+np.count_nonzero(python<=0))) if np.count_nonzero(java<=0)+np.count_nonzero(python<=0) else 1.0},
        'mawm_refinement_saved_for_followup':True}
    a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,indent=2))
    if not report['corrected_scalar_output']['exact']: raise SystemExit('Java and Python plus corrected topology fields differ')

if __name__=='__main__':main()
