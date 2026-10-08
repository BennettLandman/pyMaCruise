# Selected JIST closest-point thickness comparison

This runs the selected VandyCRUISE `ClosestPointThicknessForInner.java` mesh overload unchanged and compares its inner-surface bidirectional point data and central-surface projection with `pymacruise.thickness.closest_point_thickness` on the same saved Java `IsoSurfaceOnGrid` synthetic meshes. This holds mesh extraction fixed. The Java thickness class, `SurfaceIntersector`, `KdTree`, and `IntersectorTriangle` source bodies are compiled directly; minimal container/superclass stubs are used only for surrounding API state.

Source hash for `ClosestPointThicknessForInner.java`: `9fc19e8b9caf3defd8278a7dec1dfd417f37b0431f579da4cd17123d662534b6`.

| Field | Maximum absolute difference (mm, except vertex index) | Mean absolute difference (mm, except vertex index) |
|---|---:|---:|
| `inner.distance_inner_to_outer_mm` | 3.02761071e-07 | 5.57576937e-08 |
| `inner.distance_outer_to_inner_mm` | 0.0456767713 | 2.1310253e-05 |
| `inner.mean_thickness_mm` | 0.0228383622 | 1.06589006e-05 |
| `central.distance_inner_to_outer_mm` | 0.122684311 | 0.000204221344 |
| `central.distance_central_to_inner_mm` | 2.67322655e-07 | 3.99378497e-08 |
| `central.distance_to_inner_vertex_mm` | 0.040234823 | 5.32353038e-05 |
| `central.nearest_inner_vertex_index` | 40 | 0.0359242325 |

**Result: FAIL.** The wrapper exceeds the focused `0.001 mm` maximum-error gate: the return-distance error peaks at `0.045677 mm` and central projected-thickness error at `0.122684 mm`; 11 central nearest-inner-vertex indices differ.

At inner vertex 1050, Java selects outer face 4086 and Python face 3830. Both are equidistant (`3.610089384 mm` in Java; `3.610089337 mm` in Python), but project to opposite sides of the folded outer mesh. Returning from Java's hit gives Python `3.551857809 mm`, matching Java `3.551857696 mm`; Python's own selected hit instead leads to the reported return-distance discrepancy. This isolates the source of the largest inner error to tie ordering, not triangle-distance arithmetic.

The face-selection audit compares the chosen triangle for every inner-to-outer and central-to-inner vertex query. Face ids can differ even when projected geometry and distances agree, because adjacent triangles can share the same closest point. A highest-face-index tie rule matches Java on only about half of the mismatched queries, so the next correction must reproduce JIST tree construction and queue order rather than guess a tie policy. This comparison holds the committed Java meshes fixed and does not claim whole-workflow or real-MRI validation. See the [JSON report](jist-closest-point-thickness-java-comparison-2026-10-07.json) for counts and source hashes. Reproduce with `./.venv/bin/python scripts/compare_jist_closest_point_thickness.py`.
