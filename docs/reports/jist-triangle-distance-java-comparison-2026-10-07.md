# Selected JIST point-to-triangle distance comparison

This bounded comparison invokes VandyCRUISE's selected `KdTriangle.distance(Point3f)` method directly, scanning every target triangle for each query. Python uses the existing exact closest-point-on-triangle implementation with a centroid k-d tree. The same 41³ folded synthetic surfaces and float32 geometry are used; this isolates distance arithmetic from JIST surface extraction and the SurfaceIntersector tree traversal.

Source hashes: `KdTriangle.java` `66ae0741225bd1f57031d6122866ff5940c5975c09a98e0940a1e1810120faf1`, `KdSegment.java` `6e9c302bd828c7321a4db6ca5d9951f4c0d28b11507c63048fef40720abb8dbb`.

Queries: 256 total (64 inner vertices, 64 central vertices, 64 jittered inner vertices, 64 random interior points); target has 8136 triangles.

| Metric | Maximum error | Mean error | p99 error | Gate |
|---|---:|---:|---:|---:|
| Point distance | 4.56845854e-07 mm | 5.51019988e-08 mm | 2.31607377e-07 mm | ≤ 1e-5 mm |
| Closest-point coordinate | 0.000433490606 mm | 2.56702386e-06 mm | 2.65446557e-06 mm | ≤ 1e-3 mm |

All sampled queries pass. Java point distances match within `4.57e-7 mm`; the largest returned closest-point coordinate difference is about `0.000434 mm` (0.434 micrometres), consistent with the source's float32 geometry and arithmetic. This is direct method evidence for the point-to-triangle operation used by the thickness adapter.

The harness exhaustively scans every target triangle instead of executing `SurfaceIntersector`'s KD traversal. It does not compare JIST's mesh extractor or full thickness wrapper, and these point-query gates are specific to this synthetic check. See the machine-readable [metrics](jist-triangle-distance-java-comparison-2026-10-07.json). Reproduce with `./.venv/bin/python scripts/compare_jist_triangle_distance.py`.
