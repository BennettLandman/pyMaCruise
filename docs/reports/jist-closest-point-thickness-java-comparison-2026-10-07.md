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

**Production API result: FAIL.** It exceeds the focused `0.001 mm` maximum-error gate: return-distance error peaks at `0.045677 mm` and central projected-thickness error at `0.122684 mm`; 11 central nearest-inner-vertex indices differ.

A separate Python tree prototype now reproduces Java's complete SAH tree structure exactly on both meshes (17,707 inner-tree nodes; 32,976 outer-tree nodes). This narrows the remaining discrepancy to closest-query tie ordering: the prototype has 1 central nearest-vertex mismatch, maximum return-distance error `0.000533 mm`, and maximum projected-thickness error `0.122684 mm`. It remains an isolated prototype; production still uses the geometric centroid-tree locator.

For the largest production mismatch, at inner vertex 1050, Java selects outer face 4086 and Python face 3830. Both are equidistant (`3.610089384 mm` in Java; `3.610089337 mm` in Python), but project to opposite sides of the folded outer mesh. This comparison holds the committed Java meshes fixed and does not claim whole-workflow or real-MRI validation. See the [JSON report](jist-closest-point-thickness-java-comparison-2026-10-07.json) for source hashes and all metrics. Reproduce with `./.venv/bin/python scripts/compare_jist_closest_point_thickness.py`.
