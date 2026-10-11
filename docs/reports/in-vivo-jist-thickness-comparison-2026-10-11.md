# In-vivo JIST thickness comparison — 2026-10-11

## Scope

The selected VandyCRUISE `ClosestPointThicknessForInner` wrapper was executed on the selected Java `IsoSurfaceOnGrid` meshes generated from the real human base-workflow level sets. The Python production thickness API used those same meshes, with the same RAS-mm coordinates and Java `Point3f` float32 storage. This isolates the thickness stage; it does not compare the upstream level sets to historical MaCRUISE outputs.

## Result

**The selected wrapper comparison fails.** The main scalar fields `inner.distance_inner_to_outer_mm` and `central.distance_central_to_inner_mm` match within `2.4e-7 mm`. However, the reverse nearest-surface queries and central correspondence fields diverge at some vertices:

| Field | Maximum absolute difference | Mean absolute difference |
|---|---:|---:|
| `inner.distance_inner_to_outer_mm` | 0.000000476 mm | 0.000000068 mm |
| `inner.distance_outer_to_inner_mm` | 0.531625 mm | 0.000008682 mm |
| `inner.mean_thickness_mm` | 0.265813 mm | 0.000004356 mm |
| `central.distance_inner_to_outer_mm` | 1.695806 mm | 0.000303838 mm |
| `central.distance_central_to_inner_mm` | 0.000000231 mm | 0.000000036 mm |
| `central.distance_to_inner_vertex_mm` | 0.091162 mm | 0.000051410 mm |

The central nearest-inner-vertex index differs at **554 vertices**. The comparison script's `0.001 mm` scalar gate and exact-index gate therefore fail. These are thickness-wrapper diagnostics, not the user's surface-position acceptance metrics. Separately, the held-fixed surface extraction comparison passes the agreed Hausdorff and mean surface-position thresholds.

The human meshes contain 209 components each. For sampled discrepant queries, Java and Python return nearly equal distances but select different closest points/faces; some selected points differ by up to 2 mm. This is consistent with ambiguous nearest-surface choices among small shells, but the exact source of the remaining choice-order differences is still under investigation. No claim of thickness source equivalence is made.

The optimized Python query took **59.0 s** and peaked at **1.252 GB** process RSS. Java ran with a 2 GiB heap cap. No human images or derived image data are committed; raw geometry and outputs remain in ignored `build/`.

## Reproduction

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl .venv/bin/python \
  scripts/compare_jist_closest_point_thickness_brain_resolution.py \
  --mesh-dir build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/java-reference-surfaces/meshes \
  --levelset-dir build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/cruise-core \
  --work-dir build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/java-reference-thickness \
  --report-json build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/java-reference-thickness/human-exact-meshes-updated.json \
  --report-markdown build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/java-reference-thickness/human-exact-meshes-updated.md \
  --title 'Selected JIST closest-point thickness on in-vivo meshes' \
  --evidence-class in-vivo
```

Focused regressions after the query-kernel change: `tests/test_thickness.py`, `tests/test_jist_closest_point_thickness_report.py`, and `tests/test_jist_surface_intersector_report.py`: **12 passed in 3.23 s**. The full suite was not run in this cycle.

## Interpretation

This is a real-human-data, held-fixed stage comparison. It is neither full MaCRUISE source equivalence nor clinical validation. The next work is to reproduce Java's closest-triangle selection on the discrepant human queries and then rerun the wrapper comparison. Whole-workflow human acceptance remains open because TOADS convergence/label differences and ribbon-quality checks also remain unresolved.
