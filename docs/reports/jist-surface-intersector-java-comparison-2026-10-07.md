# Selected JIST `SurfaceIntersector` distance comparison

The harness compiles and invokes the selected JIST `SurfaceIntersector.distance(Point3f)` and its KD tree, `IntersectorTriangle`, and `KdTriangle` dependencies. Python uses `_nearest_surface_points` on the same float32 target mesh. The mesh is generated with Lewiner from the committed 41³ folded synthetic fixture; this isolates the closest-surface query from the JIST isosurface topology comparison.

Sources: `SurfaceIntersector.java` SHA-256 `0450c8cbb76770ef4458f9540980748959e0c6ec1c780f3b638173491f92916d`; `KdTree.java` `155a5e0b2a09ed1cb00a73ec99eff2288f75f6cfad55d26e659e4c72db7207f3`; `KdTriangle.java` `66ae0741225bd1f57031d6122866ff5940c5975c09a98e0940a1e1810120faf1`.

| Query coverage | Count | Max distance error (mm) | Max closest-point error (mm) |
|---|---:|---:|---:|
| Inner, central, jittered, and random interior points | 256 | 4.56845854e-07 | 0.000433490606 |

Across all 256 queries, distance error is at most `4.56845854e-07 mm`; closest-point coordinate error is at most `0.000433490606 mm`. The 1e-5 mm distance and 1e-3 mm coordinate gates pass.

This validates one selected-source nearest-surface query, not the full `ClosestPointThicknessForInner` bidirectional projection or central-surface projection. Reproduce with `./.venv/bin/python scripts/compare_jist_surface_intersector.py`; detailed metrics and source hashes are in the [JSON report](jist-surface-intersector-java-comparison-2026-10-07.json).
