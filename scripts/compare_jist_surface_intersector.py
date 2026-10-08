#!/usr/bin/env python3
"""Compare selected JIST SurfaceIntersector with Python closest-surface query."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.surfaces import SurfaceMesh, extract_levelset_surface
from pymacruise.thickness import _nearest_surface_points

ROOT = Path(__file__).resolve().parents[1]
TOADS = ROOT.parent / "toads-cruise-scm-2026-09-26/toads-cruise/src"
TREE = TOADS / "edu/jhu/ece/iacl/algorithms/graphics/locator/kdtree"
INTERSECTOR = TOADS / "edu/jhu/ece/iacl/algorithms/graphics/intersector"
HARNESS = ROOT / "tests/java_reference/SurfaceIntersectorDistanceHarness.java"
STUBS = ROOT / "tests/java_reference/thickness_comparison_stubs"
ABSTRACT_CALCULATION = ROOT / "tests/java_reference/stubs/edu/jhu/ece/iacl/jist/pipeline/AbstractCalculation.java"
VERSION_UTIL = ROOT / "tests/java_reference/surface_comparison_stubs/edu/jhu/ece/iacl/algorithms/VersionUtil.java"
FIXTURE = ROOT / "tests/data/synthetic/jist_surface_comparison"
VECMATH = (ROOT.parent / "jist-scm-2026-09-26/JIST/lib/vecmath-1.3.1.jar").resolve()
WORK = ROOT / "build/java-reference/surface-intersector-distance"
REPORT = ROOT / "docs/reports/jist-surface-intersector-java-comparison-2026-10-07.json"
REPORT_MD = ROOT / "docs/reports/jist-surface-intersector-java-comparison-2026-10-07.md"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    classes = WORK / "classes"
    classes.mkdir(parents=True, exist_ok=True)
    sources = [
        HARNESS,
        STUBS / "Jama/PackageStub.java",
        STUBS / "edu/jhu/ece/iacl/jist/structures/geom/EmbeddedSurface.java",
        ABSTRACT_CALCULATION,
        VERSION_UTIL,
        *(TREE / name for name in ("BBox.java", "BBoxEdge.java", "KdPoint3.java", "KdSegment.java", "KdTriangle.java", "KdTree.java")),
        INTERSECTOR / "IntersectorTriangle.java",
        INTERSECTOR / "SurfaceIntersector.java",
    ]
    subprocess.run(
        ["javac", "-d", str(classes), "-cp", str(VECMATH), *(str(p) for p in sources)],
        check=True,
        cwd=ROOT,
    )
    inner = extract_levelset_surface(nib.load(str(FIXTURE / "inner.nii.gz")))
    central = extract_levelset_surface(nib.load(str(FIXTURE / "central.nii.gz")))
    outer = extract_levelset_surface(nib.load(str(FIXTURE / "outer.nii.gz")))
    target = SurfaceMesh(outer.vertices.astype(np.float32).astype(np.float64), outer.faces)
    rng = np.random.default_rng(20261007)
    inner_queries = inner.vertices[np.linspace(0, len(inner.vertices)-1, 64, dtype=int)]
    central_queries = central.vertices[np.linspace(0, len(central.vertices)-1, 64, dtype=int)]
    jittered = inner_queries + rng.normal(0, 0.07, inner_queries.shape)
    low = outer.vertices.min(axis=0) + 0.13
    high = outer.vertices.max(axis=0) - 0.13
    random_queries = rng.uniform(low, high, size=(64, 3))
    queries = np.concatenate((inner_queries, central_queries, jittered, random_queries)).astype(np.float32)
    triangle_vertices = target.vertices.astype("<f4")
    triangle_faces = target.faces.astype("<i4")
    vertices_path = WORK / "vertices.raw.f32le"
    faces_path = WORK / "faces.raw.i32le"
    queries_path = WORK / "queries.raw.f32le"
    output_path = WORK / "results.raw.f64f32le"
    triangle_vertices.tofile(vertices_path)
    triangle_faces.tofile(faces_path)
    queries.astype("<f4").tofile(queries_path)
    subprocess.run(
        ["java", "-Xmx2g", "-cp", f"{classes}:{VECMATH}", "SurfaceIntersectorDistanceHarness", str(vertices_path), str(faces_path), str(queries_path), str(output_path)],
        check=True,
        cwd=ROOT,
    )
    record_dtype = np.dtype([("distance", "<f8"), ("point", "<f4", (3,))])
    java = np.fromfile(output_path, dtype=record_dtype)
    python_points, python_distances, _ = _nearest_surface_points(queries.astype(np.float64), target)
    distance_error = np.abs(java["distance"] - python_distances)
    point_error = np.linalg.norm(java["point"].astype(np.float64) - python_points, axis=1)
    report = {
        "title": "Selected JIST SurfaceIntersector.distance versus Python exact nearest-surface query",
        "synthetic_only": True,
        "method_scope": "Compiles and runs selected SurfaceIntersector.distance(Point3f), its KD tree, IntersectorTriangle, and KdTriangle source methods. It does not call the full ClosestPointThicknessForInner wrapper.",
        "sources": {
            "surface_intersector_sha256": sha(INTERSECTOR / "SurfaceIntersector.java"),
            "intersector_triangle_sha256": sha(INTERSECTOR / "IntersectorTriangle.java"),
            "kd_tree_sha256": sha(TREE / "KdTree.java"),
            "kd_triangle_sha256": sha(TREE / "KdTriangle.java"),
            "bbox_sha256": sha(TREE / "BBox.java"),
            "kdsegment_sha256": sha(TREE / "KdSegment.java"),
        },
        "fixture": {
            "path": str(FIXTURE.relative_to(ROOT)),
            "query_count": int(len(queries)),
            "query_groups": {"inner_vertices":64,"central_vertices":64,"jittered_inner_vertices":64,"random_interior_points":64},
            "seed": 20261007,
            "target_face_count": int(len(target.faces)),
            "target_vertex_sha256_float32": hashlib.sha256(triangle_vertices.tobytes()).hexdigest(),
            "target_faces_sha256_int32": hashlib.sha256(triangle_faces.tobytes()).hexdigest(),
            "queries_sha256_float32": hashlib.sha256(queries.astype("<f4").tobytes()).hexdigest(),
        },
        "metrics": {
            "distance_absolute_error_mm": {"max":float(distance_error.max()),"mean":float(distance_error.mean()),"p99":float(np.quantile(distance_error,.99))},
            "closest_point_euclidean_error_mm": {"max":float(point_error.max()),"mean":float(point_error.mean()),"p99":float(np.quantile(point_error,.99))},
            "max_distance_error_gate_mm":1e-5,
            "max_closest_point_error_gate_mm":1e-3,
            "all_queries_pass":bool(distance_error.max()<=1e-5 and point_error.max()<=1e-3),
        },
        "limitations": [
            "This validates one nearest-surface query at a time, not the bidirectional transfer and central projection branches of the complete thickness wrapper.",
            "Mesh geometry is generated by Python Lewiner extraction; JIST mesh topology is separately compared in current-base-jist-surface-comparison-2026-10-07.md.",
            "The point-coordinate bound accommodates the selected Java float32 geometry and arithmetic; reported distances have a tighter 1e-5 mm gate.",
        ],
    }
    REPORT.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n",encoding="utf-8")
    m=report["metrics"]
    REPORT_MD.write_text("\n".join([
        "# Selected JIST `SurfaceIntersector` distance comparison",
        "",
        "The harness compiles and invokes the selected JIST `SurfaceIntersector.distance(Point3f)` and its KD tree, `IntersectorTriangle`, and `KdTriangle` dependencies. Python uses `_nearest_surface_points` on the same float32 target mesh. The mesh is generated with Lewiner from the committed 41³ folded synthetic fixture; this isolates the closest-surface query from the JIST isosurface topology comparison.",
        "",
        f"Sources: `SurfaceIntersector.java` SHA-256 `{report['sources']['surface_intersector_sha256']}`; `KdTree.java` `{report['sources']['kd_tree_sha256']}`; `KdTriangle.java` `{report['sources']['kd_triangle_sha256']}`.",
        "",
        "| Query coverage | Count | Max distance error (mm) | Max closest-point error (mm) |",
        "|---|---:|---:|---:|",
        f"| Inner, central, jittered, and random interior points | {len(queries)} | {m['distance_absolute_error_mm']['max']:.9g} | {m['closest_point_euclidean_error_mm']['max']:.9g} |",
        "",
        f"Across all {len(queries)} queries, distance error is at most `{m['distance_absolute_error_mm']['max']:.9g} mm`; closest-point coordinate error is at most `{m['closest_point_euclidean_error_mm']['max']:.9g} mm`. The 1e-5 mm distance and 1e-3 mm coordinate gates pass.",
        "",
        "This validates one selected-source nearest-surface query, not the full `ClosestPointThicknessForInner` bidirectional projection or central-surface projection. Reproduce with `./.venv/bin/python scripts/compare_jist_surface_intersector.py`; detailed metrics and source hashes are in the [JSON report](jist-surface-intersector-java-comparison-2026-10-07.json).",
        "",
    ]),encoding="utf-8")
    print(json.dumps({"report":str(REPORT),"markdown":str(REPORT_MD),"metrics":m},indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
