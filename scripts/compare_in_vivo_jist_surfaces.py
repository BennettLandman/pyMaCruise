#!/usr/bin/env python3
"""Compare saved human base-workflow surfaces with selected Java JIST extraction.

This is a held-fixed level-set comparison only. It uses the three already-saved
human base level sets and the corresponding production VTK meshes. It does not
rerun the brain-resolution workflow. The Java harness runs the selected
IsoSurfaceOnGrid source on a padded crop whose boundary is strictly positive;
the pad is checked so cropping cannot remove a zero crossing.

Run from pymacruise with the project environment:

    PYTHONPATH=src .venv/bin/python scripts/compare_in_vivo_jist_surfaces.py [--core-dir PATH]

Java heap is capped at 2 GiB. Intermediate Java inputs/meshes stay in ignored
build/; the compact report is written under docs/reports/. Explicit directories
allow any saved run attempt to be compared without changing default paths.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import nibabel as nib
import numpy as np

from pymacruise.jist_intersector import JistSurfaceTree
from pymacruise.surfaces import SurfaceMesh, read_legacy_vtk


ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11"
SURFACES = CORE / "surface-thickness"
JAVA_ROOT = ROOT.parent / "toads-cruise-scm-2026-09-26/toads-cruise/src"
JIST_SOURCE = JAVA_ROOT / "edu/jhu/ece/iacl/algorithms/graphics/isosurf/IsoSurfaceOnGrid.java"
KD_SOURCE = JAVA_ROOT / "edu/jhu/ece/iacl/algorithms/graphics/locator/kdtree/KdPoint3.java"
WRAPPER = ROOT.parent / "toads-cruise-scm-2026-09-26/toads-cruise/src/edu/vanderbilt/masi/plugins/CRUISE/utilities/GetMeshSurfaceFromLevelSet.java"
JIST_LIB = ROOT.parent / "jist-scm-2026-09-26/JIST/lib"
VECMATH = JIST_LIB / "vecmath-1.3.1.jar"
STUBS = ROOT / "tests/java_reference/surface_comparison_stubs"
HARNESS = ROOT / "tests/java_reference/SurfaceIsoSurfaceOnGridHarness.java"
STUB_SOURCES = (
    STUBS / "edu/jhu/ece/iacl/jist/pipeline/AbstractCalculation.java",
    STUBS / "edu/jhu/ece/iacl/jist/structures/image/ImageData.java",
    STUBS / "edu/jhu/ece/iacl/jist/structures/image/ImageDataFloat.java",
    STUBS / "edu/jhu/ece/iacl/jist/structures/image/ImageDataDouble.java",
    STUBS / "edu/jhu/ece/iacl/jist/structures/image/ImageDataUByte.java",
    STUBS / "edu/jhu/ece/iacl/jist/structures/geom/NormalGenerator.java",
    STUBS / "edu/jhu/ece/iacl/jist/structures/geom/EmbeddedSurface.java",
    STUBS / "edu/jhu/ece/iacl/algorithms/VersionUtil.java",
    STUBS / "edu/jhu/ece/iacl/algorithms/topology/ConnectivityRule.java",
    STUBS / "javax/media/j3d/IndexedTriangleArray.java",
)
BUILD = ROOT / "build/in-vivo-jist-surface-comparison-2026-10-11"
CLASSES = BUILD / "classes"
REPORT = ROOT / "docs/reports/in-vivo-jist-surface-extraction-held-fixed-2026-10-11.md"
LEVELS = ("inner", "central", "outer")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mesh_topology(mesh: SurfaceMesh) -> dict[str, int]:
    parent = np.arange(len(mesh.vertices), dtype=np.int64)
    rank = np.zeros(len(mesh.vertices), dtype=np.uint8)
    edges: dict[tuple[int, int], int] = {}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = int(parent[x])
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        if rank[ra] == rank[rb]:
            rank[ra] += 1

    for face in mesh.faces:
        a, b, c = (int(i) for i in face)
        union(a, b)
        union(b, c)
        union(c, a)
        for u, v in ((a, b), (b, c), (c, a)):
            if u > v:
                u, v = v, u
            edges[(u, v)] = edges.get((u, v), 0) + 1
    used = np.unique(mesh.faces)
    components = len({find(int(i)) for i in used})
    return {
        "vertices": int(len(mesh.vertices)),
        "faces": int(len(mesh.faces)),
        "unique_edges": len(edges),
        "boundary_edges": sum(n == 1 for n in edges.values()),
        "nonmanifold_edges": sum(n > 2 for n in edges.values()),
        "connected_components": components,
        "euler_characteristic": int(len(mesh.vertices) - len(edges) + len(mesh.faces)),
    }


def compare_meshes(jist: SurfaceMesh, python: SurfaceMesh) -> dict[str, object]:
    # Exact selected JIST triangle query; query in bounded batches to avoid an
    # all-pairs surface-distance matrix on this brain-resolution mesh.
    def directed(source: SurfaceMesh, target: SurfaceMesh) -> np.ndarray:
        tree = JistSurfaceTree(target.vertices, target.faces)
        pieces = []
        for start in range(0, len(source.vertices), 8192):
            # query_many already returns Euclidean millimetre distances (the
            # selected Java KdTriangle.distance contract), not squared values.
            _, distances_mm, _ = tree.query_many(source.vertices[start : start + 8192])
            pieces.append(distances_mm.astype(np.float64, copy=False))
        return np.concatenate(pieces)

    d_jp = directed(jist, python)
    d_pj = directed(python, jist)
    both = np.concatenate((d_jp, d_pj))
    return {
        "jist_topology": mesh_topology(jist),
        "python_topology": mesh_topology(python),
        "symmetric_vertex_to_triangle_mm": {
            "mean": float(both.mean()),
            "hausdorff_sampled": float(both.max()),
            "jist_vertices_to_python_triangles_mean": float(d_jp.mean()),
            "python_vertices_to_jist_triangles_mean": float(d_pj.mean()),
            "jist_vertices_to_python_triangles_max": float(d_jp.max()),
            "python_vertices_to_jist_triangles_max": float(d_pj.max()),
        },
    }


def compile_java() -> None:
    required = [JIST_SOURCE, KD_SOURCE, WRAPPER, VECMATH, HARNESS, *STUB_SOURCES]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing selected JIST source or harness dependencies: " + ", ".join(missing))
    CLASSES.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["javac", "-d", str(CLASSES), "-cp", str(VECMATH), str(HARNESS),
         *(str(path) for path in STUB_SOURCES), str(KD_SOURCE), str(JIST_SOURCE)],
        cwd=ROOT, check=True,
    )


def extract_java(name: str, affine: np.ndarray) -> tuple[SurfaceMesh, dict[str, object]]:
    image_path = CORE / f"{name}.nii.gz"
    image = nib.load(str(image_path))
    values = np.asarray(image.dataobj, dtype=np.float32)
    negative = np.argwhere(values < 0.0)
    if not len(negative):
        raise ValueError(f"{name} level set has no negative samples")
    low = np.maximum(negative.min(axis=0) - 2, 0)
    high = np.minimum(negative.max(axis=0) + 3, values.shape)
    crop = np.asarray(values[tuple(slice(int(lo), int(hi)) for lo, hi in zip(low, high))], dtype="<f4")
    # The two-voxel positive halo means all crossing cubes are retained. Verify
    # the cropped boundary is strictly positive before using the smaller field.
    boundary = np.concatenate((crop[0].ravel(), crop[-1].ravel(), crop[:, 0].ravel(), crop[:, -1].ravel(), crop[:, :, 0].ravel(), crop[:, :, -1].ravel()))
    if np.any(boundary <= 0):
        raise AssertionError(f"{name} crop boundary contains nonpositive values; crop may truncate the surface")
    raw = BUILD / f"{name}.crop.f32le"
    raw.parent.mkdir(parents=True, exist_ok=True)
    crop.tofile(raw)
    out = BUILD / f"{name}.jist.json"
    heap = "-Xmx2g"
    started = time.monotonic()
    subprocess.run(
        ["java", heap, "-cp", f"{CLASSES}:{VECMATH}", "SurfaceIsoSurfaceOnGridHarness",
         str(raw), *(str(int(v)) for v in crop.shape), str(out)],
        cwd=ROOT, check=True,
    )
    runtime = time.monotonic() - started
    record = json.loads(out.read_text(encoding="utf-8"))
    if record["input_sha256"] != hashlib.sha256(crop.tobytes(order="C")).hexdigest():
        raise AssertionError(f"{name} Java extractor did not consume the exact cropped float32 field")
    vertices = np.asarray(record["vertices"], dtype=np.float64)
    vertices += low
    vertices = vertices @ affine[:3, :3].T + affine[:3, 3]
    faces = np.asarray(record["faces"], dtype=np.int64)
    # Match production's convention of positive signed volume, without altering
    # geometry/connectivity beyond winding.
    triangles = vertices[faces]
    signed = np.einsum("ij,ij->i", triangles[:, 0], np.cross(triangles[:, 1], triangles[:, 2])).sum() / 6
    if signed < 0:
        faces = faces[:, ::-1]
    info = {
        "input": str(image_path.relative_to(ROOT)),
        "input_file_sha256": sha(image_path),
        "crop_start_index": low.tolist(),
        "crop_shape": list(crop.shape),
        "crop_boundary_min": float(boundary.min()),
        "crop_boundary_max": float(boundary.max()),
        "java_input_sha256": record["input_sha256"],
        "java_runtime_seconds": runtime,
        "java_heap_limit": heap,
        "java_mesh_json_bytes": out.stat().st_size,
    }
    del values, crop, boundary, negative
    return SurfaceMesh(vertices, faces), info


def main() -> int:
    global CORE, SURFACES, BUILD, CLASSES, REPORT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-dir", type=Path, default=CORE,
                        help="directory containing inner/central/outer level-set NIfTIs")
    parser.add_argument("--surfaces-dir", type=Path, default=None,
                        help="directory containing production *_surface.vtk files (default: CORE/surface-thickness)")
    parser.add_argument("--build-dir", type=Path, default=BUILD,
                        help="ignored directory for Java classes and comparison artifacts")
    parser.add_argument("--report-path", type=Path, default=REPORT,
                        help="Markdown summary path")
    args = parser.parse_args()
    CORE = args.core_dir.resolve()
    SURFACES = (args.surfaces_dir or (CORE / "surface-thickness")).resolve()
    BUILD = args.build_dir.resolve()
    CLASSES = BUILD / "classes"
    REPORT = args.report_path.resolve()
    for name in LEVELS:
        if not (CORE / f"{name}.nii.gz").is_file():
            raise FileNotFoundError(CORE / f"{name}.nii.gz")
        if not (SURFACES / f"{name}_surface.vtk").is_file():
            raise FileNotFoundError(SURFACES / f"{name}_surface.vtk")
    BUILD.mkdir(parents=True, exist_ok=True)
    compile_java()
    report: dict[str, object] = {
        "scope": "held-fixed human level sets; direct selected Java IsoSurfaceOnGrid versus saved production Python VTK meshes",
        "core_dir": str(CORE),
        "surfaces_dir": str(SURFACES),
        "source": {
            "jist_method": "IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, CONNECT_18_6=18, iso=0, inclusive=false)",
            "isosurface_source_sha256": sha(JIST_SOURCE),
            "wrapper_source_sha256": sha(WRAPPER),
            "production_method": "PyMaCRUISE source-derived JIST connectivity-18 template extractor, saved base-surface-thickness VTK outputs",
        },
        "acceptance": {"symmetric_mean_surface_distance_mm_max": 0.01, "symmetric_sampled_hausdorff_mm_max": 1.0},
        "per_surface": {},
    }
    for name in LEVELS:
        image = nib.load(str(CORE / f"{name}.nii.gz"))
        affine = image.affine.copy()
        python_mesh = read_legacy_vtk(SURFACES / f"{name}_surface.vtk", input_frame="RAS")
        java_mesh, input_info = extract_java(name, affine)
        result = compare_meshes(java_mesh, python_mesh)
        result.update(input_info)
        result["python_vtk_sha256"] = sha(SURFACES / f"{name}_surface.vtk")
        result["python_vtk_coordinate_frame"] = python_mesh.coordinate_frame
        result["python_vtk_units"] = python_mesh.units
        result["passes_user_surface_distance_limits"] = (
            result["symmetric_vertex_to_triangle_mm"]["mean"] <= 0.01
            and result["symmetric_vertex_to_triangle_mm"]["hausdorff_sampled"] <= 1.0
        )
        report["per_surface"][name] = result
        print(json.dumps({name: result}, sort_keys=True), flush=True)

    report["all_surface_distance_limits_pass"] = all(
        result["passes_user_surface_distance_limits"]
        for result in report["per_surface"].values()
    )
    report_path = BUILD / "comparison.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [
        "# In-vivo selected JIST surface extraction comparison",
        "",
        "This compares the already-saved human base-workflow `inner`, `central`, and `outer` level sets against both the production Python VTK meshes and a direct replay of the selected Java `IsoSurfaceOnGrid.solveOriginal` source on identical float32 values. It holds the level sets fixed; it is not a full workflow run or evidence of upstream TOADS, ACE, refinement, thickness, or clinical equivalence.",
        "",
        f"Selected source SHA-256: `{report['source']['isosurface_source_sha256']}`. Wrapper SHA-256: `{report['source']['wrapper_source_sha256']}`.",
        "",
        "The Java runner uses the selected source body with small JIST container/progress/normal-generator stubs. The field is cropped around all negative samples with a two-voxel positive halo; each crop boundary is verified strictly positive, so the crop does not truncate a zero crossing. Java heap is limited to 2 GiB. Java vertices are translated by the crop start and transformed with the source NIfTI affine to RAS millimetres. The Python meshes are the saved production RAS-mm VTK outputs from the same level sets.",
        "",
        "## Geometry and topology",
        "",
        "| Level set | Java V/F/components/boundary/nonmanifold/χ | Python V/F/components/boundary/nonmanifold/χ | Symmetric vertex-to-triangle mean (mm) | Sampled Hausdorff (mm) | User limits |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for name, result in report["per_surface"].items():
        jt, pt, distance = result["jist_topology"], result["python_topology"], result["symmetric_vertex_to_triangle_mm"]
        lines.append(
            f"| {name} | {jt['vertices']}/{jt['faces']}/{jt['connected_components']}/{jt['boundary_edges']}/{jt['nonmanifold_edges']}/{jt['euler_characteristic']} | "
            f"{pt['vertices']}/{pt['faces']}/{pt['connected_components']}/{pt['boundary_edges']}/{pt['nonmanifold_edges']}/{pt['euler_characteristic']} | "
            f"{distance['mean']:.8f} | {distance['hausdorff_sampled']:.8f} | {'pass' if result['passes_user_surface_distance_limits'] else 'fail'} |"
        )
    lines.extend([
        "",
        "The reported symmetric surface statistic samples every vertex of each mesh and queries its closest point on the opposite triangle mesh using PyMaCRUISE's selected-JIST triangle-tree distance implementation. The maximum is a sampled vertex-to-triangle Hausdorff estimate, not a formal continuous triangle-to-triangle bound. Topology counts include all connected components, making small islands visible instead of silently discarding them.",
        "",
        "## Scope and limitations",
        "",
        "- This directly replays the Java `IsoSurfaceOnGrid` algorithm source, not the complete JIST runtime, VTK serializer, or historical full MaCRUISE pipeline.",
        "- Both extractors consume fixed saved human fields. Agreement validates this extractor stage only; it does not fix low GM/WM support or the fragmented output seen in the human quality audit.",
        "- Distances use the user's symmetric nearest-surface limits of 0.01 mm mean and 1 mm Hausdorff; topology agreement is reported independently and is not inferred from distance agreement.",
        f"- The machine-readable full output, hashes, crop bounds, and runtime are in `{os.path.relpath(report_path, ROOT)}` (ignored build output). The reproducible runner is [`compare_in_vivo_jist_surfaces.py`](../../scripts/compare_in_vivo_jist_surfaces.py).",
        "",
    ])
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"report": str(REPORT), "json": str(report_path), "all_surface_distance_limits_pass": report["all_surface_distance_limits_pass"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
