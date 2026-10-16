# In-vivo selected Java thickness comparison — attempt 02 — 2026-10-16

## Scope

This runs the selected VandyCRUISE `ClosestPointThicknessForInner` mesh
methods on the Java JIST meshes generated from attempt-02 human BASE level
sets. It compares those Java results with the production Python thickness
API on both the exact Java mesh coordinates and the production VTK meshes.
The selected JIST isosurface and Python meshes share the same topology and
faces; their saved coordinates differ by at most `1.34e-5 mm`, primarily from
production VTK float32 serialization. This is held-fixed mesh-stage evidence,
not full-workflow equivalence.

## Results

When Java thickness and Python thickness use the **same Java mesh
coordinates**, all seven tested outputs match bit-for-bit: inner distance in
both directions, inner mean thickness, central-to-inner distance,
central-to-inner-surface distance, central nearest-inner-vertex distance, and
nearest-inner-vertex indices.

When comparing Java thickness on Java meshes with Python thickness on the
production VTK meshes:

| Field | Mean absolute difference | Maximum absolute difference / mismatches |
| --- | ---: | ---: |
| Inner to outer distance | `2.76e-6 mm` | `0.00035 mm` |
| Outer to inner distance | `2.81e-5 mm` | `0.805 mm` |
| Inner mean thickness | `1.49e-5 mm` | `0.402 mm` |
| Central to inner distance | `2.59e-6 mm` | `0.00049 mm` |
| Central to inner surface distance | `1.43e-4 mm` | `0.1105 mm` |
| Central nearest-inner-vertex index | — | 1,962 / 240,378 differ |
| Central inner-to-outer distance | `1.397e-3 mm` | `2.282 mm` |

The discrepant thickness values are sparse. The exact same-mesh result
localizes them to the coordinate difference in the production VTK meshes and
the source nearest-surface/intersection calculation's sensitivity to tiny
coordinate changes. The production mesh extraction independently passes the
surface-position limits, but those limits do not define a thickness tolerance.
Therefore these thickness differences remain an explicit precision-sensitive
residual; they are not labeled as a pass or failure against an unprovided
tolerance.

## Reproduction

```sh
PYTHONPATH=src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_in_vivo_jist_thickness.py \
  --core-dir build/sample-singularity-base-source-profile-2026-10-15/stages/cruise-core/attempt-02 \
  --mesh-dir build/sample-singularity-base-source-profile-2026-10-15/java-surface-comparison-attempt-02 \
  --surfaces-dir build/sample-singularity-base-source-profile-2026-10-15/stages/surface-thickness/attempt-02 \
  --build-dir build/sample-singularity-base-source-profile-2026-10-15/java-thickness-comparison-attempt-02 \
  --report-path build/sample-singularity-base-source-profile-2026-10-15/java-thickness-comparison-attempt-02/comparison.json
```

The run completed under an 8 GB RSS watchdog with peak sampled RSS **1.149 GB**.
The complete metrics and source hashes are also recorded in the companion
[machine-readable report](in-vivo-root-sample-base-attempt-02-jist-thickness-2026-10-16.json).
