#!/usr/bin/env python3
"""Compare selected Java mesh thickness with production Python on human meshes.

The selected Java thickness methods are compiled from the supplied source. Both
implementations receive the same held-fixed human mesh geometry; the Python
comparison uses the production VTK meshes, whose extraction was separately
compared against JIST for the same level sets.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np

from pymacruise.jist_intersector import JistSurfaceTree
from pymacruise.surfaces import SurfaceMesh, read_legacy_vtk
from pymacruise.thickness import closest_point_thickness

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import compare_jist_closest_point_thickness as selected  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_meshes(mesh_dir: Path, surfaces_dir: Path, core_dir: Path):
    comparison = json.loads((mesh_dir / "comparison.json").read_text())
    java_meshes = {}
    python_meshes = {}
    for name in ("inner", "central", "outer"):
        payload = json.loads((mesh_dir / f"{name}.jist.json").read_text())
        # The JIST harness serializes cropped voxel coordinates. Reconstruct
        # full-volume voxel coordinates and then convert to RAS millimetres,
        # exactly as the matched surface comparator does.
        vertices = np.asarray(payload["vertices"], dtype=np.float64)
        vertices += np.asarray(comparison["per_surface"][name]["crop_start_index"])
        affine = nib.load(str(core_dir / f"{name}.nii.gz")).affine
        vertices = vertices @ affine[:3, :3].T + affine[:3, 3]
        java_meshes[name] = SurfaceMesh(
            vertices.astype(np.float32).astype(np.float64),
            np.asarray(payload["faces"], dtype=np.int64),
            coordinate_frame="RAS",
            units="mm",
        )
        python_meshes[name] = read_legacy_vtk(
            surfaces_dir / f"{name}_surface.vtk", input_frame="RAS"
        )
        if len(java_meshes[name].vertices) != len(python_meshes[name].vertices):
            raise ValueError(f"{name} Java/Python vertex counts differ")
        if not np.array_equal(java_meshes[name].faces, python_meshes[name].faces):
            raise ValueError(f"{name} Java/Python ordered faces differ")
    return java_meshes, python_meshes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-dir", type=Path, required=True,
                        help="directory containing inner/central/outer level-set NIfTIs")
    parser.add_argument("--mesh-dir", type=Path, required=True,
                        help="directory containing inner/central/outer.jist.json")
    parser.add_argument("--surfaces-dir", type=Path, required=True,
                        help="directory containing production *_surface.vtk files")
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--report-path", type=Path, required=True)
    args = parser.parse_args()
    core_dir = args.core_dir.resolve()
    mesh_dir, surfaces_dir = args.mesh_dir.resolve(), args.surfaces_dir.resolve()
    build_dir, report_path = args.build_dir.resolve(), args.report_path.resolve()
    classes = build_dir / "classes"
    classes.mkdir(parents=True, exist_ok=True)
    toads = selected.TOADS
    tree = selected.TREE
    intersector = selected.INTERSECTOR
    sources = [
        selected.HARNESS, selected.DIAGNOSTIC, selected.FACE_HARNESS,
        selected.TREE_HARNESS, selected.STUB / "Jama/PackageStub.java",
        selected.STUB / "edu/jhu/ece/iacl/jist/structures/geom/EmbeddedSurface.java",
        selected.STUB / "edu/jhu/ece/iacl/algorithms/thickness/ThicknessSolver.java",
        selected.STUB / "edu/jhu/ece/iacl/jist/structures/image/ImageData.java",
        selected.STUB / "edu/jhu/ece/iacl/algorithms/graphics/isosurf/IsoSurfaceOnGrid.java",
        selected.STUB / "edu/jhu/ece/iacl/algorithms/topology/ConnectivityRule.java",
        selected.ABSTRACT, selected.VERSION,
        *[tree / name for name in ("BBox.java", "BBoxEdge.java", "KdPoint3.java",
                                  "KdSegment.java", "KdTriangle.java", "KdTree.java")],
        intersector / "IntersectorTriangle.java", intersector / "SurfaceIntersector.java",
        selected.SELECTED,
    ]
    subprocess.run(["javac", "-d", str(classes), "-cp", str(selected.VECMATH),
                    *map(str, sources)], check=True, cwd=ROOT)

    java_meshes, python_meshes = load_meshes(mesh_dir, surfaces_dir, core_dir)
    raw_paths = {}
    for name, mesh in java_meshes.items():
        vertices = build_dir / f"{name}.vertices.f32le"
        faces = build_dir / f"{name}.faces.i32le"
        mesh.vertices.astype("<f4").tofile(vertices)
        mesh.faces.astype("<i4").tofile(faces)
        raw_paths[name] = (vertices, faces)
    result_path = build_dir / "thickness-java.raw"
    subprocess.run([
        "java", "-Xmx2g", "-cp", f"{classes}:{selected.VECMATH}",
        "ClosestPointThicknessForInnerHarness",
        str(raw_paths["inner"][0]), str(raw_paths["inner"][1]),
        str(raw_paths["central"][0]), str(raw_paths["central"][1]),
        str(raw_paths["outer"][0]), str(raw_paths["outer"][1]), str(result_path),
    ], check=True, cwd=ROOT)

    java_raw = result_path.read_bytes()
    offset = 0
    inner_count = int(np.frombuffer(java_raw, dtype="<i4", count=1, offset=offset)[0]); offset += 4
    java_inner = np.frombuffer(java_raw, dtype="<f8", count=inner_count * 3,
                               offset=offset).reshape(inner_count, 3).copy()
    offset += inner_count * 24
    central_count = int(np.frombuffer(java_raw, dtype="<i4", count=1, offset=offset)[0]); offset += 4
    java_central = np.frombuffer(java_raw, dtype="<f8", count=central_count * 4,
                                 offset=offset).reshape(central_count, 4).copy()
    if inner_count != len(java_meshes["inner"].vertices) or central_count != len(java_meshes["central"].vertices):
        raise ValueError("Java thickness output counts do not match input meshes")

    python = closest_point_thickness(
        python_meshes["inner"], python_meshes["outer"], central=python_meshes["central"]
    )
    python_on_java_meshes = closest_point_thickness(
        java_meshes["inner"], java_meshes["outer"], central=java_meshes["central"]
    )
    vertex_geometry_max_mm = {
        name: float(np.max(np.linalg.norm(
            java_meshes[name].vertices - python_meshes[name].vertices, axis=1
        )))
        for name in java_meshes
    }
    inner_names = ("distance_inner_to_outer_mm", "distance_outer_to_inner_mm", "mean_thickness_mm")
    central_names = ("distance_inner_to_outer_mm", "distance_central_to_inner_mm",
                     "distance_to_inner_vertex_mm", "nearest_inner_vertex_index")
    metrics = {}
    same_mesh_metrics = {}
    for index, name in enumerate(inner_names):
        expected = python.inner_surface.point_data[name]
        error = np.abs(java_inner[:, index] - expected)
        metrics[f"inner.{name}"] = {"max_abs_difference": float(error.max()),
                                      "mean_abs_difference": float(error.mean())}
        same_error = np.abs(
            java_inner[:, index] - python_on_java_meshes.inner_surface.point_data[name]
        )
        same_mesh_metrics[f"inner.{name}"] = {
            "max_abs_difference": float(same_error.max()),
            "mean_abs_difference": float(same_error.mean()),
        }
    for index, name in enumerate(central_names):
        expected = python.central_surface.point_data[name]
        if name == "nearest_inner_vertex_index":
            metrics[f"central.{name}"] = {
                "mismatches": int(np.count_nonzero(java_central[:, index] != expected)),
                "query_count": int(len(expected)),
            }
        else:
            error = np.abs(java_central[:, index] - expected)
            metrics[f"central.{name}"] = {"max_abs_difference": float(error.max()),
                                            "mean_abs_difference": float(error.mean())}
            same_error = np.abs(
                java_central[:, index]
                - python_on_java_meshes.central_surface.point_data[name]
            )
            same_mesh_metrics[f"central.{name}"] = {
                "max_abs_difference": float(same_error.max()),
                "mean_abs_difference": float(same_error.mean()),
            }
    same_mesh_metrics["central.nearest_inner_vertex_index"] = {
        "mismatches": int(np.count_nonzero(
            java_central[:, 3]
            != python_on_java_meshes.central_surface.point_data["nearest_inner_vertex_index"]
        )),
        "query_count": int(central_count),
    }
    report = {
        "scope": "selected Java ClosestPointThicknessForInner source vs production Python on held-fixed human meshes",
        "synthetic_only": False,
        "mesh_dir": str(mesh_dir),
        "surfaces_dir": str(surfaces_dir),
        "core_dir": str(core_dir),
        "source_sha256": {
            "ClosestPointThicknessForInner": sha(selected.SELECTED),
            "SurfaceIntersector": sha(selected.INTERSECTOR / "SurfaceIntersector.java"),
            "IntersectorTriangle": sha(selected.INTERSECTOR / "IntersectorTriangle.java"),
            "KdTree": sha(selected.TREE / "KdTree.java"),
            "KdTriangle": sha(selected.TREE / "KdTriangle.java"),
        },
        "mesh_counts": {name: {"vertices": len(mesh.vertices), "faces": len(mesh.faces)}
                        for name, mesh in java_meshes.items()},
        "corresponding_vertex_geometry_max_mm": vertex_geometry_max_mm,
        "metrics": metrics,
        "same_java_mesh_python_metrics": same_mesh_metrics,
        "python_summary": {"mean_thickness_mm": python.mean_mm,
                           "standard_deviation_mm": python.standard_deviation_mm,
                           "nonzero_inner_vertices": python.nonzero_vertex_count},
        "python_on_java_mesh_summary": {
            "mean_thickness_mm": python_on_java_meshes.mean_mm,
            "standard_deviation_mm": python_on_java_meshes.standard_deviation_mm,
            "nonzero_inner_vertices": python_on_java_meshes.nonzero_vertex_count,
        },
        "java_heap_limit": "-Xmx2g",
        "limitations": [
            "This is held-fixed human mesh thickness evidence, not whole-workflow equivalence.",
            "Python production VTK geometry was independently compared with Java JIST on the same human level sets.",
            "The source wrapper's serialized VTK/NIfTI outputs are not exercised.",
        ],
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"report": str(report_path), "metrics": metrics,
                      "python_summary": report["python_summary"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
