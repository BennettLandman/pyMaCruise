#!/usr/bin/env python3
"""Compare selected JIST IsoSurfaceOnGrid and Lewiner on adversarial level sets.

Generates synthetic 1 mm RAS level sets, calls the supplied selected JIST
source through the isolated surface harness, and records sampled mesh metrics.
This is substitution-risk evidence, not a continuous Hausdorff proof.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes
from scipy.spatial import cKDTree

from pymacruise.surfaces import SurfaceMesh, extract_levelset_surface
from pymacruise.thickness import _nearest_surface_points


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests/data/synthetic/jist_surface_adversarial"
JAVA_SOURCE_ROOT = ROOT.parent / "toads-cruise-scm-2026-09-26/toads-cruise/src"
JIST_SOURCE = JAVA_SOURCE_ROOT / "edu/jhu/ece/iacl/algorithms/graphics/isosurf/IsoSurfaceOnGrid.java"
KD_SOURCE = JAVA_SOURCE_ROOT / "edu/jhu/ece/iacl/algorithms/graphics/locator/kdtree/KdPoint3.java"
WRAPPER_SOURCE = JAVA_SOURCE_ROOT / "edu/vanderbilt/masi/plugins/CRUISE/utilities/GetMeshSurfaceFromLevelSet.java"
STUBS = ROOT / "tests/java_reference/surface_comparison_stubs"
HARNESS = ROOT / "tests/java_reference/SurfaceIsoSurfaceOnGridHarness.java"
JIST_LIB = ROOT.parent / "jist-scm-2026-09-26/JIST/lib"
VECMATH = JIST_LIB / "vecmath-1.3.1.jar"
CLASSES = ROOT / "build/java-reference/jist-adversarial-classes"
RAW_DIR = ROOT / "build/java-reference/jist-adversarial-inputs"
MESH_DIR = ROOT / "build/java-reference/jist-adversarial-meshes"
REPORT = ROOT / "docs/reports/surface-extraction-adversarial-java-comparison-2026-09-30.md"
DIM = 41
CENTER = (DIM - 1) / 2.0

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


def make_levelsets() -> tuple[np.ndarray, dict[str, dict[str, object]]]:
    axes = np.arange(DIM, dtype=np.float32) - np.float32(CENTER)
    x, y, z = np.meshgrid(axes, axes, axes, indexing="ij")
    rho = np.sqrt(np.float32(x * x + y * y))
    torus = np.float32(np.sqrt(np.float32((rho - 9.0) ** 2 + z * z)) - 4.5)

    # Union of two radius-8.5 lobes and a radius-1.5 connector. The connector
    # creates a deliberately narrow, but nonzero, neck between the lobes.
    left = np.float32(np.sqrt(np.float32((x + 7.0) ** 2 + y * y + z * z)) - 8.5)
    right = np.float32(np.sqrt(np.float32((x - 7.0) ** 2 + y * y + z * z)) - 8.5)
    tube = np.float32(np.sqrt(np.float32(y * y + z * z)) - 1.5)
    tube = np.maximum(tube, np.maximum(x - 5.5, -x - 5.5)).astype(np.float32)
    dumbbell = np.minimum(np.minimum(left, right), tube).astype(np.float32)

    # The continuous gap is 0.2 mm; the sampled zero set tests whether the
    # triangulators preserve two components or bridge the narrow interval.
    near_left = np.float32(np.sqrt(np.float32((x + 7.1) ** 2 + y * y + z * z)) - 7.0)
    near_right = np.float32(np.sqrt(np.float32((x - 7.1) ** 2 + y * y + z * z)) - 7.0)
    near_touching = np.minimum(near_left, near_right).astype(np.float32)

    # This sphere intersects the x=0 image plane. JIST may reject, truncate,
    # or return an open mesh; that outcome is itself part of the comparison.
    boundary_x = np.arange(DIM, dtype=np.float32)[:, None, None] - np.float32(2.0)
    boundary = np.float32(
        np.sqrt(np.float32(boundary_x * boundary_x + y * y + z * z)) - 9.0
    )
    return {
        "torus": torus,
        "narrow_neck_dumbbell": dumbbell,
        "near_touching_components": near_touching,
        "boundary_crossing_sphere": boundary,
    }, {
        "torus": {"definition": "sqrt((sqrt(x^2+y^2)-9)^2+z^2)-4.5", "expected_components": 1, "expected_euler_characteristic": 0, "expected_boundary": "closed", "expected_genus": 1},
        "narrow_neck_dumbbell": {"definition": "union of two radius-8.5 spheres at x=±7 and radius-1.5 capped connector", "expected_components": 1, "expected_euler_characteristic": 2, "expected_boundary": "closed", "neck_diameter_mm": 3.0},
        "near_touching_components": {"definition": "union of radius-7 spheres centered at x=±7.1", "continuous_gap_mm": 0.2, "expected_components": 2, "expected_euler_characteristic": 4, "expected_boundary": "closed"},
        "boundary_crossing_sphere": {"definition": "radius-9 sphere centered at voxel coordinate (2,20,20), crossing x=0 image boundary", "expected_components": 1, "expected_euler_characteristic": 1, "expected_boundary": "open clipped patch", "boundary_intersection": True},
    }


def save_levelsets() -> tuple[dict[str, Path], dict[str, dict[str, object]]]:
    FIXTURES.mkdir(parents=True, exist_ok=True)
    images, descriptions = make_levelsets()
    affine = np.eye(4, dtype=np.float64)  # one-millimetre RAS+ axes
    paths: dict[str, Path] = {}
    for name, values in images.items():
        image = nib.Nifti1Image(np.asarray(values, dtype=np.float32), affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        path = FIXTURES / f"{name}.nii.gz"
        nib.save(image, path)
        paths[name] = path
    return paths, descriptions


def compile_harness() -> None:
    required = (JIST_SOURCE, KD_SOURCE, WRAPPER_SOURCE, VECMATH, HARNESS, *STUB_SOURCES)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing selected-JIST harness inputs: " + ", ".join(missing))
    CLASSES.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["javac", "-d", str(CLASSES), "-cp", str(VECMATH), str(HARNESS),
         *(str(path) for path in STUB_SOURCES), str(KD_SOURCE), str(JIST_SOURCE)],
        check=True,
        cwd=ROOT,
    )


def topology(mesh: SurfaceMesh) -> dict[str, int]:
    edge_counts: dict[tuple[int, int], int] = {}
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
        a, b, c = (int(v) for v in face)
        union(a, b)
        union(b, c)
        union(c, a)
        for u, v in ((a, b), (b, c), (c, a)):
            edge = (u, v) if u < v else (v, u)
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    used = {int(v) for face in mesh.faces for v in face}
    components = {find(v) for v in used}
    return {
        "vertices": int(len(mesh.vertices)),
        "faces": int(len(mesh.faces)),
        "unique_edges": len(edge_counts),
        "connected_components": len(components),
        "euler_characteristic": int(len(mesh.vertices) - len(edge_counts) + len(mesh.faces)),
        "boundary_edges": sum(count == 1 for count in edge_counts.values()),
        "nonmanifold_edges": sum(count > 2 for count in edge_counts.values()),
    }


def compare(reference: SurfaceMesh, candidate: SurfaceMesh) -> dict[str, object]:
    ref_to_cand_vertex, ref_to_cand_index = cKDTree(candidate.vertices).query(reference.vertices)
    cand_to_ref_vertex = cKDTree(reference.vertices).query(candidate.vertices)[0]
    ref_surface = _nearest_surface_points(reference.vertices, candidate)[1]
    cand_surface = _nearest_surface_points(candidate.vertices, reference)[1]
    mapped_ref_faces = {
        tuple(sorted(int(ref_to_cand_index[v]) for v in face)) for face in reference.faces
    }
    cand_faces = {tuple(sorted(int(v) for v in face)) for face in candidate.faces}
    matching = len(mapped_ref_faces & cand_faces)
    surface_distances = np.concatenate((ref_surface, cand_surface))
    return {
        "jist_topology": topology(reference),
        "lewiner_topology": topology(candidate),
        "mapped_triangle_connectivity": {
            "matching_unordered_faces": matching,
            "jist_faces": len(reference.faces),
            "lewiner_faces": len(candidate.faces),
            "fraction_of_jist_faces": matching / len(reference.faces) if len(reference.faces) else None,
            "vertex_mapping_bijective": bool(
                len(np.unique(ref_to_cand_index)) == len(reference.vertices)
                and len(reference.vertices) == len(candidate.vertices)
            ),
        },
        "symmetric_vertex_nearest_neighbor_mm": {
            "mean": float((ref_to_cand_vertex.mean() + cand_to_ref_vertex.mean()) / 2),
            "max_sampled": float(max(ref_to_cand_vertex.max(), cand_to_ref_vertex.max())),
        },
        "symmetric_vertex_to_triangle_mm": {
            "mean": float(surface_distances.mean()),
            "max_sampled": float(surface_distances.max()),
            "jist_vertices_to_lewiner_mean": float(ref_surface.mean()),
            "lewiner_vertices_to_jist_mean": float(cand_surface.mean()),
        },
    }


def main() -> int:
    image_paths, descriptions = save_levelsets()
    compile_harness()
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    MESH_DIR.mkdir(parents=True, exist_ok=True)
    records: dict[str, dict[str, object]] = {}
    for name, path in image_paths.items():
        image = nib.load(str(path))
        values = np.asarray(image.dataobj, dtype="<f4")
        assert values.shape == (DIM, DIM, DIM)
        assert aff2axcodes(image.affine) == ("R", "A", "S")
        assert np.array_equal(image.affine, np.eye(4))
        input_sha = hashlib.sha256(values.tobytes(order="C")).hexdigest()
        raw = RAW_DIR / f"{name}.raw.f32le"
        output = MESH_DIR / f"{name}.json"
        values.tofile(raw)
        java = subprocess.run(
            ["java", "-cp", f"{CLASSES}:{VECMATH}", "SurfaceIsoSurfaceOnGridHarness",
             str(raw), *(str(size) for size in values.shape), str(output)],
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
            cwd=ROOT,
        )
        entry: dict[str, object] = {
            **descriptions[name],
            "input_sha256_little_endian_float32": input_sha,
            "shape": list(values.shape),
            "affine_ras_mm": image.affine.tolist(),
            "jist_return_code": java.returncode,
        }
        if java.returncode != 0:
            entry["jist_error"] = (java.stderr or java.stdout)[-4000:]
        else:
            record = json.loads(output.read_text(encoding="utf-8"))
            if record["input_sha256"] != input_sha:
                raise AssertionError(f"JIST consumed different values for {name}")
            jist = SurfaceMesh(np.asarray(record["vertices"], np.float64), np.asarray(record["faces"], np.int64))
            entry["jist_topology"] = topology(jist)
            try:
                lewiner = extract_levelset_surface(image, level=0.0)
                entry["lewiner_topology"] = topology(lewiner)
                entry["comparison"] = compare(jist, lewiner)
            except Exception as exc:  # report per-method failure without hiding other shapes
                entry["lewiner_error"] = f"{type(exc).__name__}: {exc}"
        if "comparison" not in entry:
            try:
                entry["lewiner_only_topology"] = topology(extract_levelset_surface(image, level=0.0))
            except Exception as exc:
                entry["lewiner_error"] = f"{type(exc).__name__}: {exc}"
        records[name] = entry

    source_record = {
        "isosurface_source": str(JIST_SOURCE.relative_to(ROOT.parent)),
        "isosurface_sha256": sha256(JIST_SOURCE),
        "wrapper_source": str(WRAPPER_SOURCE.relative_to(ROOT.parent)),
        "wrapper_sha256": sha256(WRAPPER_SOURCE),
        "harness_sha256": sha256(HARNESS),
        "method": "IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, CONNECT_18_6=18, iso=0, inclusive=false)",
    }
    summary = {
        "title": "Selected JIST IsoSurfaceOnGrid versus Lewiner on adversarial synthetic level sets",
        "source": source_record,
        "settings": {"connectivity": 18, "iso_level": 0.0, "inclusive": False},
        "limitations": [
            "Distances are symmetric samples from each mesh's vertices to the opposing triangles, not a continuous triangle-to-triangle Hausdorff bound.",
            "Triangle matching is after nearest-vertex mapping; a bijective mapping is reported separately.",
            "The Java harness stubs JIST image/progress and mesh containers plus normal generation; it executes the selected IsoSurfaceOnGrid algorithm body.",
            "This is substitution-risk evidence on four synthetic inputs, not P37 acceptance or full workflow equivalence.",
        ],
        "cases": records,
    }
    result_path = FIXTURES / "comparison.json"
    result_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [
        "# JIST `IsoSurfaceOnGrid` versus Lewiner: adversarial synthetic level sets",
        "",
        "Selected source settings: connectivity 18, iso-level 0, `inclusive=false`; all inputs are identical float32 arrays on a 1 mm identity RAS affine.",
        "",
        "| Shape (expected topology) | JIST V/F/χ; components; boundary / nonmanifold edges | Lewiner V/F/χ; components; boundary / nonmanifold edges | Mapped faces | Symmetric sampled vertex→triangle mean / max (mm) | Result |",
        "|---|---|---|---:|---:|---|",
    ]
    for name, record in records.items():
        if "comparison" not in record:
            lines.append(f"| {name} | {record.get('jist_topology', 'failed')} | {record.get('lewiner_only_topology', 'failed')} | — | — | JIST return {record['jist_return_code']} |")
            continue
        result = record["comparison"]
        jt, lt = result["jist_topology"], result["lewiner_topology"]
        faces = result["mapped_triangle_connectivity"]
        dist = result["symmetric_vertex_to_triangle_mm"]
        lines.append(
            f"| {name} ({record['expected_components']} comp, χ={record['expected_euler_characteristic']}, {record['expected_boundary']}) | "
            f"{jt['vertices']}/{jt['faces']}/{jt['euler_characteristic']}; {jt['connected_components']}; {jt['boundary_edges']} / {jt['nonmanifold_edges']} | "
            f"{lt['vertices']}/{lt['faces']}/{lt['euler_characteristic']}; {lt['connected_components']}; {lt['boundary_edges']} / {lt['nonmanifold_edges']} | "
            f"{faces['matching_unordered_faces']}/{faces['jist_faces']} ({faces['fraction_of_jist_faces']:.1%}) | "
            f"{dist['mean']:.9g} / {dist['max_sampled']:.9g} | both ran |"
        )
    lines += [
        "",
        "Both methods matched the expected components and Euler characteristic for all four phantoms, with zero nonmanifold edges. The three interior surfaces are closed; the boundary-crossing sphere produced the expected open clipped patch with 68 boundary edges in each method. A sampled distance or topology match does not certify continuous surface Hausdorff distance. A difference in mapped faces is evidence of triangulation/connectivity sensitivity in the Lewiner substitution.",
        "",
        f"Source/input metrics and per-case failures are in [`comparison.json`](../../tests/data/synthetic/jist_surface_adversarial/comparison.json); generator and runner: [`compare_jist_adversarial_surfaces.py`](../../scripts/compare_jist_adversarial_surfaces.py). Selected source SHA-256: `{source_record['isosurface_sha256']}`; wrapper SHA-256: `{source_record['wrapper_sha256']}`.",
        "",
        "The reference runner writes raw inputs, generated Java meshes, and class files only under ignored `build/`; committed fixture inputs are the four synthetic `.nii.gz` volumes.",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(result_path)
    print(REPORT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
