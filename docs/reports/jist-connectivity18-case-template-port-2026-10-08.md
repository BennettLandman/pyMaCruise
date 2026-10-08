# JIST connectivity-18 surface templates: source comparison (2026-10-08)

The production extractor now uses a Python-native 256-case table derived by executing the supplied JIST `IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, 18, 0, false)` source. It preserves the selected source’s zero/near-zero preprocessing, case connectivity, edge interpolation, shared grid-edge vertices, face index ordering, and float32 output points. No Java runtime or original source checkout is needed at runtime. The table is packaged at [`jist_isosurface_connectivity18_cases.json`](../../src/pymacruise/data/jist_isosurface_connectivity18_cases.json).

## Evidence

- Selected source: `IsoSurfaceOnGrid.java`, CVS r1.7; SHA-256 `a323ac5eed8b66ad62a68cd93a6bc7bcf01281453fc6541ee5c4ccc529d40c38`.
- Development-only Java case harness: `tests/java_reference/SurfaceIsoSurfaceOnGridCaseMatrixHarness.java`; SHA-256 `b81ce8ab5de4eee42981db06fda8a3116863411f4aafd0dcbb68c14a1934d6da`.
- Every one of 256 sign patterns was generated using five distinct deterministic magnitude assignments. The emitted edge vertices and triangle templates were invariant across all five scalar assignments: 1,536 edge-vertex occurrences and 820 triangles across the five-case matrices.
- Focused pytest command: `./.venv/bin/pytest -q tests/test_jist_connectivity18_templates.py tests/test_surfaces.py tests/test_plus_thickness.py tests/test_plus_file_workflow.py tests/test_jist_triangle_distance_report.py tests/test_workflow_preflight_cli.py` — **278 passed, 0 failed, 6 expected RAS-orientation warnings in 26.81s**.
- Workflow regression command: `./.venv/bin/pytest -q tests/test_base_workflow_restartable.py tests/test_plus_file_workflow.py tests/test_plus_full_workflow.py` — **8 passed in 607.55s**, peak sampled RSS **0.638 GB** under a 16 GB watchdog. This reruns base restart and plus full/file workflow integration against the new extractor.
- The Java comparison harness consumed the exact saved float32 level sets for the synthetic base and plus workflows. Both comparisons used an 8 GB RSS watchdog; peak sampled RSS was 0.535 GB (base) and 0.128 GB (plus).
- Packaging check: built a wheel and installed it into `build/jist18-wheel-check/site`; from outside the checkout, the installed module found the packaged JSON and extracted a synthetic sphere (192 vertices, 380 faces, finite RAS-mm coordinates). Runtime dependencies no longer include scikit-image; it remains a development dependency for independent reference fixtures.

## Full workflow level-set comparisons

| Saved synthetic workflow | Level set | JIST V/F | Python V/F | Components / Euler (both) | Ordered faces matched | Paired vertex mean / max (mm) |
|---|---|---:|---:|---:|---:|---:|
| base | inner | 116166/232328 | 116166/232328 | 1 / 2 | 232328/232328 | 1.24e-06 / 7.14e-06 |
| base | central | 121074/242144 | 121074/242144 | 1 / 2 | 242144/242144 | 1.26e-06 / 7.14e-06 |
| base | outer | 123828/247652 | 123828/247652 | 1 / 2 | 247652/247652 | 1.23e-06 / 7.14e-06 |
| plus | inner | 3724/7444 | 3724/7444 | 1 / 2 | 7444/7444 | 1.39e-07 / 8.72e-07 |
| plus | central | 3724/7444 | 3724/7444 | 1 / 2 | 7444/7444 | 1.38e-07 / 8.72e-07 |
| plus | outer | 3724/7444 | 3724/7444 | 1 / 2 | 7444/7444 | 1.36e-07 / 8.72e-07 |

All six saved level sets have bijective nearest-vertex mappings and exact ordered triangle correspondence after mapping. Since corresponding triangles have the same connectivity, the maximum paired-vertex displacement bounds the continuous distance between each pair of planar triangles; the largest bound across the six meshes is **7.15×10⁻⁶ mm**, below the agreed 0.01 mm mean-position tolerance and 1 mm Hausdorff tolerance. All meshes are closed, manifold, RAS-mm, and match JIST components and Euler characteristics exactly.

## Scope

This is held-fixed selected-source equivalence for the surface extraction stage on the saved synthetic base and plus level sets, plus exhaustive local case-template validation. It does not establish that upstream workflow outputs match MATLAB/JIST, validate the available human MRI, or establish paper-level correctness of the selected historical JIST implementation. The old Lewiner comparison reports remain historical evidence and are superseded for current production surfaces.

Commands used for the saved-workflow comparisons:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 --interval-seconds 0.1 --report-interval-seconds 30 -- ./.venv/bin/python scripts/compare_workflow_jist_surfaces.py --core-dir build/ace-indexed-heap-2026-10-08/base-seven-stage-run/stages/cruise-core/attempt-01 --work-dir build/base-seven-stage-jist18-source-comparison-2026-10-08
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 --interval-seconds 0.1 --report-interval-seconds 30 -- ./.venv/bin/python scripts/compare_workflow_jist_surfaces.py --core-dir build/wheel-plus-workflow-current/run-wheel-validated/plus-surface/stages/surfaces_thickness/attempt-01 --work-dir build/plus-seven-stage-jist18-source-comparison-2026-10-08
```
