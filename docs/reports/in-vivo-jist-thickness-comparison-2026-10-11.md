# In-vivo JIST thickness comparison — 2026-10-11

## Scope

The selected VandyCRUISE `ClosestPointThicknessForInner` wrapper was executed on the selected Java `IsoSurfaceOnGrid` meshes generated from the real human base-workflow level sets. The Python production thickness API used those same meshes, with the same RAS-mm coordinates and Java `Point3f` float32 storage. This isolates the thickness stage; it does not compare the upstream level sets to historical MaCRUISE outputs.

## Result

**The selected wrapper comparison passes exactly after a source-precision fix.** The previous comparison had up to 1.695806 mm central-field error and 554 nearest-inner-vertex index mismatches. Inspection showed that Java's `Math.sqrt(float)` promotes its input to double, while Numba had selected a float32 `math.sqrt` overload because the input remained float32. Explicitly promoting the triangle squared-distance value to float64 before `sqrt` restored Java's double-precision square root and its nearest-face ordering. The box-distance path now mirrors `javax.vecmath.Point3f.distance`: float32 square/sum, double square root, then float32 narrowing.

| Field | Maximum absolute difference | Mean absolute difference |
|---|---:|---:|
| `inner.distance_inner_to_outer_mm` | 0 mm | 0 mm |
| `inner.distance_outer_to_inner_mm` | 0 mm | 0 mm |
| `inner.mean_thickness_mm` | 0 mm | 0 mm |
| `central.distance_inner_to_outer_mm` | 0 mm | 0 mm |
| `central.distance_central_to_inner_mm` | 0 mm | 0 mm |
| `central.distance_to_inner_vertex_mm` | 0 mm | 0 mm |
| `central.nearest_inner_vertex_index` | 0 | 0 |

The central nearest-inner-vertex index now matches at every vertex (**0 mismatches**); all seven fields pass the comparison script's `0.001 mm` scalar gate and exact-index gate. These are thickness-wrapper diagnostics, separate from the user's surface-position acceptance metrics. The held-fixed surface extraction comparison independently passes the agreed Hausdorff and mean surface-position thresholds.

The held-fixed meshes contain 209 components each; this thickness comparison does not establish their anatomical validity or resolve the upstream ribbon-quality deficit. Mesh extraction and upstream level sets remain separate validation questions.

The corrected Python query took **61.2 s** and peaked at **1.257 GB** process RSS. Java ran with a 2 GiB heap cap. No human images or derived image data are committed; raw geometry and outputs remain in ignored `build/`.

## Reproduction

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl .venv/bin/python \
  scripts/compare_jist_closest_point_thickness_brain_resolution.py \
  --mesh-dir build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/java-reference-surfaces/meshes \
  --levelset-dir build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/cruise-core \
  --work-dir build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/java-reference-thickness \
  --report-json build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/java-reference-thickness/human-double-sqrt.json \
  --report-markdown build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/java-reference-thickness/human-double-sqrt.md \
  --title 'Selected JIST closest-point thickness on in-vivo meshes' \
  --evidence-class in-vivo
```

Focused regressions after the precision correction, plus the saved Java TGDM source-match regression, ran on code committed as `6b51612`: `tests/test_thickness.py`, `tests/test_jist_closest_point_thickness_report.py`, `tests/test_jist_surface_intersector_report.py`, and `tests/test_base_core_current_java_match.py`: **15 passed in 56.31 s**. The full suite was not run in this cycle.

## Interpretation

This is a real-human-data, held-fixed source-stage comparison. It demonstrates exact agreement for the selected thickness wrapper on the same Java meshes. It does not establish full MaCRUISE source equivalence or clinical validity. Whole-workflow human acceptance remains open because TOADS convergence/label differences and ribbon-quality checks remain unresolved.
