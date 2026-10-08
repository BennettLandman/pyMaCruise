# Selected JIST closest-point thickness comparison

This runs the selected VandyCRUISE `ClosestPointThicknessForInner.java` mesh overload unchanged and compares its inner-surface bidirectional point data and central-surface projection with `pymacruise.thickness.closest_point_thickness` on the same float32-stored synthetic meshes. The Java class, `SurfaceIntersector`, `KdTree`, and `IntersectorTriangle` source bodies are compiled directly; minimal container/superclass stubs are used only for surrounding API state.

Source hash for `ClosestPointThicknessForInner.java`: `9fc19e8b9caf3defd8278a7dec1dfd417f37b0431f579da4cd17123d662534b6`.

| Field | Maximum absolute difference (mm, except vertex index) | Mean absolute difference (mm, except vertex index) |
|---|---:|---:|
| `inner.distance_inner_to_outer_mm` | 2.4818026e-07 | 5.58960313e-08 |
| `inner.distance_outer_to_inner_mm` | 0.0456766371 | 2.12665298e-05 |
| `inner.mean_thickness_mm` | 0.0228382951 | 1.06371544e-05 |
| `central.distance_inner_to_outer_mm` | 0.307235328 | 0.000459474508 |
| `central.distance_central_to_inner_mm` | 2.13801647e-07 | 3.70309591e-08 |
| `central.distance_to_inner_vertex_mm` | 0.0500429311 | 0.000127196019 |
| `central.nearest_inner_vertex_index` | 134 | 0.343239713 |

**Result: FAIL.** The wrapper exceeds the focused `0.001 mm` maximum-error gate: the return-distance error peaks at `0.045677 mm` and central projected-thickness error at `0.307235 mm`; 21 central nearest-inner-vertex indices differ.

At inner vertex 1124, Java selects outer face 3717 and Python face 3714. Both are equidistant (`3.610089384 mm` in Java; `3.610089337 mm` in Python), but project to opposite sides of the folded outer mesh. Returning from Java's hit gives Python `3.551857809 mm`, matching Java `3.551857830 mm`; Python's own selected hit instead leads to the reported return-distance discrepancy. This isolates the source of the largest inner error to tie ordering, not triangle-distance arithmetic.

The current production operator remains a Python geometric implementation and does not reproduce Java's float32 triangle arithmetic and KD traversal ordering at exact ties. The report does not establish workflow equivalence or real-MRI behavior. It uses Lewiner meshes to isolate thickness behavior; surface topology remains a separate P37 deviation. The compact output with source hashes and mesh counts is in the [JSON report](jist-closest-point-thickness-java-comparison-2026-10-07.json). Reproduce with `./.venv/bin/python scripts/compare_jist_closest_point_thickness.py`.
