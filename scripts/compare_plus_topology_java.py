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
    def read_java(name):
        with gzip.open(root/name,'rb') as f:
            header=tuple(int.from_bytes(f.read(4),'big') for _ in range(3))
            values=np.frombuffer(f.read(),dtype='>f4').astype(np.float32)
        if header!=shape or values.size!=int(np.prod(shape)): raise SystemExit(f'{name}: Java output dimensions/payload mismatch')
        return values.reshape(shape,order='F')
    topology_java=read_java('topology-java.raw.gz')
    topology_python=np.load(root/'python-corrected.npy').astype(np.float32)
    topo_diff=np.abs(topology_java.astype(np.float64)-topology_python.astype(np.float64))
    eroded_java=read_java('topology-mawm-eroded.raw.gz')
    eroded_python=np.load(root/'python-mawm-eroded.npy').astype(np.float32)
    erode_diff=np.abs(eroded_java.astype(np.float64)-eroded_python.astype(np.float64))
    refined_java=read_java('topology-mawm-refined.raw.gz')
    refined_python=np.load(root/'python-atlas_refined.npy').astype(np.float32)
    refined_diff=np.abs(refined_java.astype(np.float64)-refined_python.astype(np.float64))
    sampled_rss=float(m.group(1))
    records={
        'topology_corrected':{'exact':bool(np.array_equal(topology_java,topology_python)),'changed_voxels':int(np.count_nonzero(topology_java!=topology_python)),'voxel_count':int(topology_java.size),'maximum_absolute_error':float(topo_diff.max(initial=0)),'mean_absolute_error':float(topo_diff.mean())},
        'mawm_eroded':{'exact':bool(np.array_equal(eroded_java,eroded_python)),'changed_voxels':int(np.count_nonzero(eroded_java!=eroded_python)),'voxel_count':int(eroded_java.size),'maximum_absolute_error':float(erode_diff.max(initial=0)),'mean_absolute_error':float(erode_diff.mean())},
        'mawm_refined':{'exact':bool(np.array_equal(refined_java,refined_python)),'changed_voxels':int(np.count_nonzero(refined_java!=refined_python)),'voxel_count':int(refined_java.size),'maximum_absolute_error':float(refined_diff.max(initial=0)),'mean_absolute_error':float(refined_diff.mean())},
    }
    report={
        'synthetic_only':True,'comparison':'selected Java AlgorithmTopologyCorrection + MASImorph refinement vs plus Python topology_correct','scope':provenance['scope'],
        'workflow_run_signature':provenance['workflow_run_signature'],'shape_xyz':list(shape),'coordinate_system':'RAS','java_peak_sampled_rss_gb':sampled_rss if sampled_rss>=0.01 else None,'java_rss_limit_gb':8,
        'rss_measurement_note':None if sampled_rss>=0.01 else 'The watchdog returned less than 0.01 GB for this short JVM run; treated as below measurement resolution, not a reliable peak.',
        'outputs':records,'all_outputs_exact':all(row['exact'] for row in records.values())}
    a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,indent=2))
    if not report['all_outputs_exact']: raise SystemExit('Java and Python plus topology/refinement outputs differ')

if __name__=='__main__':main()
