# Selected JIST closest-point thickness comparison

This runs the selected VandyCRUISE `ClosestPointThicknessForInner.java` mesh overload unchanged and compares its inner-surface bidirectional point data and central-surface projection with `pymacruise.thickness.closest_point_thickness` on the same saved Java `IsoSurfaceOnGrid` synthetic meshes. This holds mesh extraction fixed. The Java thickness class, `SurfaceIntersector`, `KdTree`, and `IntersectorTriangle` source bodies are compiled directly; minimal container/superclass stubs are used only for surrounding API state.

Source hash for `ClosestPointThicknessForInner.java`: `9fc19e8b9caf3defd8278a7dec1dfd417f37b0431f579da4cd17123d662534b6`.

| Field | Maximum absolute difference (mm, except vertex index) | Mean absolute difference (mm, except vertex index) |
|---|---:|---:|
| `inner.distance_inner_to_outer_mm` | 0 | 0 |
| `inner.distance_outer_to_inner_mm` | 0 | 0 |
| `inner.mean_thickness_mm` | 0 | 0 |
| `central.distance_inner_to_outer_mm` | 0 | 0 |
| `central.distance_central_to_inner_mm` | 0 | 0 |
| `central.distance_to_inner_vertex_mm` | 2.03233304e-07 | 3.88426888e-08 |
| `central.nearest_inner_vertex_index` | 0 | 0 |

**Production API result: PASS.** On these held-fixed Java meshes, central nearest-vertex indices match exactly, five of six physical fields are bit-exact, and the remaining field differs by less than `2.1e-7 mm`; physical-field maximum error is `2.03233304e-07 mm` and central vertex-index mismatches are 0.

The connected Python JIST tree reproduces Java's complete SAH tree structure exactly on both meshes (17,707 inner-tree nodes; 32,976 outer-tree nodes). Its face-query selections and wrapper output metrics are exact on this fixture. The older centroid-tree locator is retained only as a baseline: it selects different faces at ties (569 inner-to-outer and 932 central-to-inner query mismatches).

This is selected-source synthetic mesh evidence only. It does not validate upstream processing, mesh extraction equivalence, whole-workflow outputs, or real MRI. Exact JIST tree construction increases runtime; brain-resolution cost remains unmeasured. See the [JSON report](jist-closest-point-thickness-java-comparison-2026-10-07.json) for source hashes and metrics. Reproduce with `./.venv/bin/python scripts/compare_jist_closest_point_thickness.py`.
