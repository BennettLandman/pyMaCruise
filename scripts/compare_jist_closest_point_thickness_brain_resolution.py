#!/usr/bin/env python3
"""Compare selected JIST thickness wrapper on the committed full-grid meshes.

First generate the Java IsoSurfaceOnGrid meshes with
``compare_jist_surface_current_base.py``. Mesh arrays are kept under ignored
``build/``; only hashes and comparison metrics are committed.
"""
from __future__ import annotations

import hashlib
import argparse
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import nibabel as nib
import numpy as np

from pymacruise.surfaces import SurfaceMesh
from pymacruise.thickness import closest_point_thickness
import compare_jist_closest_point_thickness as selected


ROOT = Path(__file__).resolve().parents[1]
MESHES = ROOT / "build/java-reference/current-base-surface-comparison/meshes"
LEVELSETS = ROOT / "tests/data/synthetic/base_core_source_match/java_reference"
REPORT = ROOT / "docs/reports/jist-closest-point-thickness-brain-resolution-java-comparison-2026-10-07.json"
REPORT_MD = ROOT / "docs/reports/jist-closest-point-thickness-brain-resolution-java-comparison-2026-10-07.md"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-dir", type=Path, default=MESHES,
                        help="directory containing inner/central/outer Java mesh JSON")
    parser.add_argument("--levelset-dir", type=Path, default=LEVELSETS,
                        help="directory containing matching inner/central/outer NIfTIs")
    parser.add_argument("--work-dir", type=Path, default=selected.WORK,
                        help="ignored build directory for compiled harness and arrays")
    parser.add_argument("--report-json", type=Path, default=REPORT)
    parser.add_argument("--report-markdown", type=Path, default=REPORT_MD)
    parser.add_argument("--title", default="Selected JIST closest-point thickness on synthetic meshes")
    parser.add_argument(
        "--evidence-class", choices=("synthetic", "in-vivo"), default="synthetic",
        help="label the source-comparison dataset correctly in the JSON and Markdown reports",
    )
    args = parser.parse_args()
    mesh_dir = args.mesh_dir.resolve()
    levelset_dir = args.levelset_dir.resolve()
    work_dir = args.work_dir.resolve()
    report_json = args.report_json.resolve()
    report_markdown = args.report_markdown.resolve()
    classes = work_dir / "classes"
    classes.mkdir(parents=True, exist_ok=True)
    sources = [
        selected.HARNESS,
        selected.STUB / "Jama/PackageStub.java",
        selected.STUB / "edu/jhu/ece/iacl/jist/structures/geom/EmbeddedSurface.java",
        selected.STUB / "edu/jhu/ece/iacl/algorithms/thickness/ThicknessSolver.java",
        selected.STUB / "edu/jhu/ece/iacl/jist/structures/image/ImageData.java",
        selected.STUB / "edu/jhu/ece/iacl/algorithms/graphics/isosurf/IsoSurfaceOnGrid.java",
        selected.STUB / "edu/jhu/ece/iacl/algorithms/topology/ConnectivityRule.java",
        selected.ABSTRACT,
        selected.VERSION,
        *[
            selected.TREE / name
            for name in ("BBox.java", "BBoxEdge.java", "KdPoint3.java", "KdSegment.java", "KdTriangle.java", "KdTree.java")
        ],
        selected.INTERSECTOR / "IntersectorTriangle.java",
        selected.INTERSECTOR / "SurfaceIntersector.java",
        selected.SELECTED,
    ]
    subprocess.run(
        ["javac", "-d", str(classes), "-cp", str(selected.VECMATH), *map(str, sources)],
        check=True,
        cwd=ROOT,
    )
    mesh_paths = {name: mesh_dir / f"{name}.json" for name in ("inner", "central", "outer")}
    missing = [str(path) for path in mesh_paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Java meshes are missing; first run "
            "./.venv/bin/python scripts/compare_workflow_jist_surfaces.py"
        )

    meshes = {}
    for name, path in mesh_paths.items():
        payload = json.loads(path.read_text(encoding="utf-8"))
        index_vertices = np.asarray(payload["vertices"], dtype=np.float32).astype(np.float64)
        image = nib.load(str(levelset_dir / f"{name}.nii.gz"))
        vertices_ras = nib.affines.apply_affine(image.affine, index_vertices)
        meshes[name] = SurfaceMesh(vertices_ras, np.asarray(payload["faces"], dtype=np.int64))

    raw_paths = {}
    for name, mesh in meshes.items():
        vertices_path = work_dir / f"{name}-vertices.raw.f32le"
        faces_path = work_dir / f"{name}-faces.raw.i32le"
        mesh.vertices.astype("<f4").tofile(vertices_path)
        mesh.faces.astype("<i4").tofile(faces_path)
        raw_paths[name] = (vertices_path, faces_path)

    java_output = work_dir / "result.raw"
    subprocess.run(
        [
            "java", "-Xmx2g", "-cp", f"{classes}:{selected.VECMATH}",
            "ClosestPointThicknessForInnerHarness",
            str(raw_paths["inner"][0]), str(raw_paths["inner"][1]),
            str(raw_paths["central"][0]), str(raw_paths["central"][1]),
            str(raw_paths["outer"][0]), str(raw_paths["outer"][1]), str(java_output),
        ],
        check=True,
        cwd=ROOT,
    )

    started = time.perf_counter()
    python_result = closest_point_thickness(
        meshes["inner"], meshes["outer"], central=meshes["central"]
    )
    python_seconds = time.perf_counter() - started
    raw = java_output.read_bytes()
    offset = 0
    inner_count = int(np.frombuffer(raw, dtype="<i4", count=1, offset=offset)[0])
    offset += 4
    java_inner = np.frombuffer(raw, dtype="<f8", count=inner_count * 3, offset=offset).reshape(inner_count, 3).copy()
    offset += inner_count * 24
    central_count = int(np.frombuffer(raw, dtype="<i4", count=1, offset=offset)[0])
    offset += 4
    java_central = np.frombuffer(raw, dtype="<f8", count=central_count * 4, offset=offset).reshape(central_count, 4).copy()

    errors = {}
    inner_names = ("distance_inner_to_outer_mm", "distance_outer_to_inner_mm", "mean_thickness_mm")
    central_names = ("distance_inner_to_outer_mm", "distance_central_to_inner_mm", "distance_to_inner_vertex_mm", "nearest_inner_vertex_index")
    for index, name in enumerate(inner_names):
        errors[f"inner.{name}"] = np.abs(java_inner[:, index] - python_result.inner_surface.point_data[name])
    for index, name in enumerate(central_names):
        errors[f"central.{name}"] = np.abs(java_central[:, index] - python_result.central_surface.point_data[name])
    maximum = {name: float(values.max()) for name, values in errors.items()}
    mean = {name: float(values.mean()) for name, values in errors.items()}
    physical_error = max(value for name, value in maximum.items() if name.endswith("_mm"))
    index_mismatches = int(
        np.count_nonzero(
            java_central[:, 3]
            != python_result.central_surface.point_data["nearest_inner_vertex_index"]
        )
    )
    peak_rss_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (
        1e9 if sys.platform == "darwin" else 1e6
    )
    result = {
        "title": args.title,
        "synthetic_only": args.evidence_class == "synthetic",
        "human_data": args.evidence_class == "in-vivo",
        "method_scope": "Compiles the selected VandyCRUISE ClosestPointThicknessForInner mesh overload and selected JIST tree/intersector sources unchanged. One full-grid Java wrapper execution is compared with the production Python API on identical Java IsoSurfaceOnGrid meshes.",
        "sources": {
            "closest_point_thickness_for_inner_sha256": sha(selected.SELECTED),
            "surface_intersector_sha256": sha(selected.INTERSECTOR / "SurfaceIntersector.java"),
            "intersector_triangle_sha256": sha(selected.INTERSECTOR / "IntersectorTriangle.java"),
            "kd_tree_sha256": sha(selected.TREE / "KdTree.java"),
            "kd_triangle_sha256": sha(selected.TREE / "KdTriangle.java"),
        },
        "fixture": {
            "levelset_path": str(levelset_dir),
            "mesh_path": str(mesh_dir),
            "mesh_origin": "selected Java IsoSurfaceOnGrid outputs; mesh extraction is held fixed",
            "levelset_sha256": {name: sha(levelset_dir / f"{name}.nii.gz") for name in meshes},
            "mesh_json_sha256": {name: sha(path) for name, path in mesh_paths.items()},
            "mesh_counts": {name: {"vertices": len(mesh.vertices), "faces": len(mesh.faces)} for name, mesh in meshes.items()},
            "coordinate_storage": "RAS millimetres stored as Java Point3f float32",
        },
        "metrics": {
            "max_absolute_field_difference": maximum,
            "mean_absolute_field_difference": mean,
            "maximum_physical_field_error_gate_mm": 0.001,
            "central_nearest_inner_vertex_index_mismatches": index_mismatches,
            "all_physical_fields_pass": physical_error <= 0.001,
            "nearest_vertex_indices_exact": index_mismatches == 0,
            "comparison_pass": physical_error <= 0.001 and index_mismatches == 0,
        },
        "runtime": {
            "python_api_elapsed_seconds": python_seconds,
            "python_process_peak_rss_gb": peak_rss_gb,
            "java_heap_cap_gb": 2,
        },
        "limitations": [
            "This validates the selected thickness wrapper on held-fixed Java meshes, not the complete workflow.",
            "Mesh extraction connectivity and upstream stage equivalence remain separate P37 items.",
            "Python peak RSS is measured for this process; Java heap is capped at 2 GiB, but process RSS is not included.",
        ],
    }
    report_json.parent.mkdir(parents=True, exist_ok=True)
    report_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [
        f"# {args.title}",
        "",
        "The selected VandyCRUISE thickness wrapper and production Python API ran on the same Java `IsoSurfaceOnGrid` meshes. This isolates the wrapper; upstream processing and mesh extraction are held fixed.",
        "",
        "| Mesh | Vertices | Faces |",
        "|---|---:|---:|",
        *[f"| {name} | {len(mesh.vertices):,} | {len(mesh.faces):,} |" for name, mesh in meshes.items()],
        "",
        "| Field | Maximum absolute difference | Mean absolute difference |",
        "|---|---:|---:|",
        *[f"| `{name}` | {maximum[name]:.9g} | {mean[name]:.9g} |" for name in maximum],
        "",
        f"**Selected wrapper comparison: {'PASS' if result['metrics']['comparison_pass'] else 'FAIL'}.** Maximum physical-field error is `{physical_error:.9g} mm` against the `0.001 mm` scalar gate; central nearest-vertex index mismatches: {index_mismatches}.",
        "",
        f"The Python API took `{python_seconds:.1f} s` and peaked at `{peak_rss_gb:.3f} GB` process RSS. Java ran with a 2 GiB heap cap. This is selected-source {args.evidence_class} stage evidence, not whole-workflow equivalence. Mesh extraction connectivity differences are addressed separately. See the [JSON report]({report_json.name}) for source and fixture hashes. Reproduce with `./.venv/bin/python scripts/compare_jist_closest_point_thickness_brain_resolution.py --mesh-dir {mesh_dir} --levelset-dir {levelset_dir} --work-dir {work_dir} --report-json {report_json} --report-markdown {report_markdown} --evidence-class {args.evidence_class}`.",
        "",
    ]
    report_markdown.parent.mkdir(parents=True, exist_ok=True)
    report_markdown.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"report": str(report_json), "metrics": result["metrics"], "runtime": result["runtime"]}, indent=2))
    return 0 if result["metrics"]["comparison_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
