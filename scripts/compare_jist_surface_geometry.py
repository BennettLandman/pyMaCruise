#!/usr/bin/env python3
"""Compare selected Java IsoSurfaceOnGrid with the Python Lewiner extractor.

Run with the project environment and the source checkouts in their documented
sibling locations:

    .venv/bin/python scripts/compare_jist_surface_geometry.py

The report and Java mesh arrays are written under tests/data/synthetic and
docs/reports; Java build classes and raw temporary inputs stay under ignored
build/.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np
from scipy.spatial import cKDTree

from pymacruise.surfaces import SurfaceMesh, extract_levelset_surface
from pymacruise.thickness import _nearest_surface_points, closest_point_thickness


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "data" / "synthetic" / "jist_surface_comparison"
JAVA_SOURCE_ROOT = ROOT.parent / "toads-cruise-scm-2026-09-26" / "toads-cruise" / "src"
JIST_SOURCE = JAVA_SOURCE_ROOT / "edu/jhu/ece/iacl/algorithms/graphics/isosurf/IsoSurfaceOnGrid.java"
KD_SOURCE = JAVA_SOURCE_ROOT / "edu/jhu/ece/iacl/algorithms/graphics/locator/kdtree/KdPoint3.java"
JIST_WRAPPER_SOURCE = ROOT.parent / "toads-cruise-scm-2026-09-26" / "toads-cruise" / "src/edu/vanderbilt/masi/plugins/CRUISE/utilities/GetMeshSurfaceFromLevelSet.java"
STUBS = ROOT / "tests" / "java_reference" / "surface_comparison_stubs"
HARNESS = ROOT / "tests" / "java_reference" / "SurfaceIsoSurfaceOnGridHarness.java"
CLASSES = ROOT / "build" / "java-reference" / "surface-comparison-classes"
REPORT = ROOT / "docs" / "reports" / "surface-extraction-java-comparison-2026-09-30.md"
JAVA_MESH_DIR = FIXTURES / "java_meshes"
RAW_DIR = ROOT / "build" / "java-reference" / "surface-comparison-inputs"

JIST_LIB = ROOT.parent / "jist-scm-2026-09-26" / "JIST" / "lib"
VECMATH = JIST_LIB / "vecmath-1.3.1.jar"

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


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_java() -> dict[str, Path]:
    required = (JIST_SOURCE, KD_SOURCE, JIST_WRAPPER_SOURCE, VECMATH, HARNESS)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Selected JIST source/reference runtime is unavailable; missing: "
            + ", ".join(missing)
        )
    if not all(path.is_file() for path in STUB_SOURCES):
        raise FileNotFoundError("The direct Java harness is missing one or more API stubs")

    CLASSES.mkdir(parents=True, exist_ok=True)
    JAVA_MESH_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    compile_command = [
        "javac", "-d", str(CLASSES), "-cp", str(VECMATH),
        str(HARNESS), *(str(path) for path in STUB_SOURCES),
        str(KD_SOURCE), str(JIST_SOURCE),
    ]
    subprocess.run(compile_command, check=True, cwd=ROOT)

    result: dict[str, Path] = {}
    for name in ("inner", "central", "outer"):
        image_path = FIXTURES / f"{name}.nii.gz"
        image = nib.load(str(image_path))
        values = np.asarray(image.dataobj, dtype="<f4")
        raw_path = RAW_DIR / f"{name}.raw.f32le"
        values.tofile(raw_path)
        json_path = JAVA_MESH_DIR / f"{name}.json"
        subprocess.run(
            [
                "java", "-cp", f"{CLASSES}:{VECMATH}",
                "SurfaceIsoSurfaceOnGridHarness", str(raw_path),
                *(str(size) for size in values.shape), str(json_path),
            ],
            check=True,
            cwd=ROOT,
        )
        result[name] = json_path
    return result


def mesh_topology(mesh: SurfaceMesh) -> dict[str, int]:
    edges: dict[tuple[int, int], int] = {}
    parent = list(range(len(mesh.vertices)))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for face in mesh.faces:
        a, b, c = (int(value) for value in face)
        union(a, b)
        union(b, c)
        union(c, a)
        for u, v in ((a, b), (b, c), (c, a)):
            key = (u, v) if u < v else (v, u)
            edges[key] = edges.get(key, 0) + 1

    used = {int(vertex) for face in mesh.faces for vertex in face}
    components = {find(value) for value in used}
    return {
        "vertices": int(len(mesh.vertices)),
        "faces": int(len(mesh.faces)),
        "unique_edges": int(len(edges)),
        "boundary_edges": int(sum(count == 1 for count in edges.values())),
        "nonmanifold_edges": int(sum(count > 2 for count in edges.values())),
        "connected_components": int(len(components)),
        "euler_characteristic": int(len(mesh.vertices) - len(edges) + len(mesh.faces)),
    }


def compare_pair(reference: SurfaceMesh, candidate: SurfaceMesh) -> dict[str, object]:
    ref_to_candidate_vertex, ref_to_candidate_index = cKDTree(candidate.vertices).query(
        reference.vertices
    )
    candidate_to_ref_vertex = cKDTree(reference.vertices).query(candidate.vertices)[0]
    ref_to_candidate_surface = _nearest_surface_points(reference.vertices, candidate)[1]
    candidate_to_ref_surface = _nearest_surface_points(candidate.vertices, reference)[1]
    symmetric_surface = np.concatenate((ref_to_candidate_surface, candidate_to_ref_surface))
    mapped_reference_faces = {
        tuple(sorted(int(ref_to_candidate_index[vertex]) for vertex in face))
        for face in reference.faces
    }
    candidate_faces = {
        tuple(sorted(int(vertex) for vertex in face)) for face in candidate.faces
    }
    matching_faces = len(mapped_reference_faces & candidate_faces)
    return {
        "jist_topology": mesh_topology(reference),
        "lewiner_topology": mesh_topology(candidate),
        "triangle_connectivity_after_nearest_vertex_mapping": {
            "matching_unordered_faces": matching_faces,
            "jist_face_count": len(reference.faces),
            "lewiner_face_count": len(candidate.faces),
            "matching_fraction_of_jist_faces": float(matching_faces / len(reference.faces)),
            "vertex_mapping_is_bijective": bool(
                len(np.unique(ref_to_candidate_index)) == len(reference.vertices)
                and len(reference.vertices) == len(candidate.vertices)
            ),
        },
        "vertex_nearest_neighbor_mm": {
            "symmetric_mean": float((ref_to_candidate_vertex.mean() + candidate_to_ref_vertex.mean()) / 2),
            "symmetric_max": float(max(ref_to_candidate_vertex.max(), candidate_to_ref_vertex.max())),
            "jist_vertices_to_lewiner_mean": float(ref_to_candidate_vertex.mean()),
            "lewiner_vertices_to_jist_mean": float(candidate_to_ref_vertex.mean()),
        },
        "vertex_to_triangle_surface_distance_mm": {
            "symmetric_mean": float(symmetric_surface.mean()),
            "symmetric_max": float(symmetric_surface.max()),
            "jist_vertices_to_lewiner_mean": float(ref_to_candidate_surface.mean()),
            "lewiner_vertices_to_jist_mean": float(candidate_to_ref_surface.mean()),
        },
    }


def load_java_mesh(path: Path) -> tuple[dict[str, object], SurfaceMesh]:
    record = json.loads(path.read_text(encoding="utf-8"))
    mesh = SurfaceMesh(
        np.asarray(record["vertices"], dtype=np.float64),
        np.asarray(record["faces"], dtype=np.int64),
    )
    return record, mesh


def main() -> int:
    java_mesh_paths = run_java()
    comparisons: dict[str, dict[str, object]] = {}
    jist_meshes: dict[str, SurfaceMesh] = {}
    lewiner_meshes: dict[str, SurfaceMesh] = {}
    source_input_hashes: dict[str, str] = {}

    for name, java_path in java_mesh_paths.items():
        image_path = FIXTURES / f"{name}.nii.gz"
        image = nib.load(str(image_path))
        values = np.asarray(image.dataobj, dtype="<f4")
        record, jist_mesh = load_java_mesh(java_path)
        actual_hash = hashlib.sha256(values.tobytes(order="C")).hexdigest()
        if record["input_sha256"] != actual_hash:
            raise AssertionError(f"Java harness did not consume the exact {name} float32 input")
        if tuple(record["dimensions"]) != tuple(values.shape):
            raise AssertionError(f"Java harness dimension mismatch for {name}")
        source_input_hashes[name] = actual_hash
        lewiner_mesh = extract_levelset_surface(image, level=0.0)
        jist_meshes[name] = jist_mesh
        lewiner_meshes[name] = lewiner_mesh
        comparisons[name] = compare_pair(jist_mesh, lewiner_mesh)

    jist_thickness = closest_point_thickness(
        jist_meshes["inner"], jist_meshes["outer"], central=jist_meshes["central"]
    )
    lewiner_thickness = closest_point_thickness(
        lewiner_meshes["inner"], lewiner_meshes["outer"], central=lewiner_meshes["central"]
    )
    thickness = {
        "measurement": "same PyMaCRUISE closest_point_thickness operator on each triangulator's nested meshes",
        "jist_geometry_mean_mm": float(jist_thickness.mean_mm),
        "lewiner_geometry_mean_mm": float(lewiner_thickness.mean_mm),
        "mean_delta_mm": float(lewiner_thickness.mean_mm - jist_thickness.mean_mm),
        "jist_geometry_sd_mm": float(jist_thickness.standard_deviation_mm),
        "lewiner_geometry_sd_mm": float(lewiner_thickness.standard_deviation_mm),
        "jist_nonzero_inner_vertices": int(jist_thickness.nonzero_vertex_count),
        "lewiner_nonzero_inner_vertices": int(lewiner_thickness.nonzero_vertex_count),
    }
    summary = {
        "title": "Selected JIST IsoSurfaceOnGrid versus scikit-image Lewiner",
        "source": {
            "isosurface_source": str(JIST_SOURCE.relative_to(ROOT.parent)),
            "isosurface_sha256": sha256(JIST_SOURCE),
            "wrapper_source": str(JIST_WRAPPER_SOURCE.relative_to(ROOT.parent)),
            "wrapper_sha256": sha256(JIST_WRAPPER_SOURCE),
            "method": "IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, CONNECT_18_6=18, iso=0, inclusive=false)",
            "version": "IsoSurfaceOnGrid CVS r1.7; GetMeshSurfaceFromLevelSet CVS r1.4",
            "input_hashes_little_endian_float32": source_input_hashes,
        },
        "fixture": "41x41x41 1 mm RAS cortex-like folded sphere; inner/central/outer signed level sets; see tests/data/synthetic/jist_surface_comparison/README.md",
        "coordinate_note": "Both mesh outputs are in voxel-index coordinates for this identity-affine fixture; these are numerically millimetres.",
        "per_surface": comparisons,
        "thickness": thickness,
        "limitations": [
            "The Java harness compiles the selected IsoSurfaceOnGrid source directly but stubs JIST progress/image containers and the normal generator; it compares returned point/triangle geometry, not the complete JIST runtime or VTK writer.",
            "Vertex-nearest distances depend on tessellation density. The directed vertex-to-triangle distances sample each mesh's vertices against the other mesh's triangles; they are not a continuous triangle-to-triangle Hausdorff bound.",
            "The vertex coordinates and topology counts align closely, but only the reported fraction of unordered faces match after nearest-vertex mapping; differing local connectivity can affect interpolated triangle interiors and mesh-dependent measurements.",
            "Thickness uses the same Python closest-point operator on both sets of meshes; this isolates geometry impact but does not compare Java thickness code.",
            "One folded synthetic fixture does not establish clinical, general, or end-to-end equivalence.",
        ],
    }
    result_path = FIXTURES / "comparison.json"
    result_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_lines = [
        "# JIST `IsoSurfaceOnGrid` versus Lewiner: folded synthetic phantom",
        "",
        "This report compares the selected Java `IsoSurfaceOnGrid.solveOriginal` implementation with PyMaCRUISE's scikit-image Lewiner extraction on the exact same three nested float32 level sets. The source wrapper selects connectivity 18, iso-level 0, and `inclusive=false`.",
        "",
        f"Source `IsoSurfaceOnGrid.java` SHA-256: `{summary['source']['isosurface_sha256']}`. Wrapper SHA-256: `{summary['source']['wrapper_sha256']}`.",
        "",
        "## Surface metrics",
        "",
        "| Level set | JIST V/F/χ | Lewiner V/F/χ | JIST boundary/nonmanifold edges | Lewiner boundary/nonmanifold edges | Symmetric vertex NN mean / max (mm) | Symmetric vertex-to-triangle mean / max (mm) |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for name, result in comparisons.items():
        jt = result["jist_topology"]
        lt = result["lewiner_topology"]
        vnn = result["vertex_nearest_neighbor_mm"]
        vts = result["vertex_to_triangle_surface_distance_mm"]
        report_lines.append(
            f"| {name} | {jt['vertices']}/{jt['faces']}/{jt['euler_characteristic']} | "
            f"{lt['vertices']}/{lt['faces']}/{lt['euler_characteristic']} | "
            f"{jt['boundary_edges']}/{jt['nonmanifold_edges']} | "
            f"{lt['boundary_edges']}/{lt['nonmanifold_edges']} | "
            f"{vnn['symmetric_mean']:.6f} / {vnn['symmetric_max']:.6f} | "
            f"{vts['symmetric_mean']:.6f} / {vts['symmetric_max']:.6f} |"
        )
    report_lines.extend(
        [
            "",
            "The vertex-count/topology invariants match on this phantom, but this does not mean the triangulations match. After bijectively mapping near-identical Java vertices to Lewiner vertex indices, unordered triangles match by the fractions below. The remaining face connectivity can change interpolated triangle interiors and mesh-dependent measurements.",
            "",
            "| Level set | Matching unordered triangles | Matching fraction |",
            "|---|---:|---:|",
        ]
    )
    for name, result in comparisons.items():
        connectivity = result["triangle_connectivity_after_nearest_vertex_mapping"]
        report_lines.append(
            f"| {name} | {connectivity['matching_unordered_faces']}/{connectivity['jist_face_count']} | "
            f"{connectivity['matching_fraction_of_jist_faces']:.1%} |"
        )
    report_lines.extend(
        [
            "",
            "## Downstream thickness impact",
            "",
            f"Using the same `closest_point_thickness` operator on the two nested mesh sets: JIST geometry mean thickness `{thickness['jist_geometry_mean_mm']:.6f} mm` (SD `{thickness['jist_geometry_sd_mm']:.6f}`, {thickness['jist_nonzero_inner_vertices']} nonzero inner vertices); Lewiner geometry mean thickness `{thickness['lewiner_geometry_mean_mm']:.6f} mm` (SD `{thickness['lewiner_geometry_sd_mm']:.6f}`, {thickness['lewiner_nonzero_inner_vertices']} nonzero inner vertices); difference `{thickness['mean_delta_mm']:+.6f} mm`.",
            "",
            "## Scope and limitations",
            "",
            "- The harness compiles and invokes the selected `IsoSurfaceOnGrid` source directly. Minimal API stubs replace JIST progress/image containers and normal generation; the harness records the Java algorithm's actual returned vertices and triangle indices.",
            "- Both methods consume identical float32 values. Per-volume SHA-256 digests are checked by the runner. The fixture's identity affine makes index coordinates equal to RAS millimetres.",
            "- Vertex-nearest distances depend on tessellation density. Vertex-to-triangle distances here sample mesh vertices against the other mesh's triangles; they are not a continuous triangle-to-triangle Hausdorff bound.",
            "- Although vertex coordinates and topology counts align closely, only the reported fraction of unordered faces match after mapping vertices; this indicates differing local connectivity/interpolation and is a relevant substitution risk for mesh-dependent measurements.",
            "- Thickness uses the same Python closest-point operator on both extracted mesh sets. This quantifies geometry's effect on that measure; it is not a direct Java thickness implementation comparison.",
            "- This single synthetic folded fixture does not establish clinical or end-to-end equivalence. Matching this case would not remove the need for broader shapes, adversarial topology cases, and source-version verification.",
            "",
            "Full metric output and source/input provenance: [`comparison.json`](../../tests/data/synthetic/jist_surface_comparison/comparison.json). Fixture generator: [`make_phantom.py`](../../tests/data/synthetic/jist_surface_comparison/make_phantom.py).",
        ]
    )
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(json.dumps({"comparison": str(result_path), "report": str(REPORT), "per_surface": comparisons, "thickness": thickness}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
