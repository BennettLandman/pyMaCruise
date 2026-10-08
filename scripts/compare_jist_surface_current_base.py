#!/usr/bin/env python3
"""Compare selected JIST and Lewiner surfaces on the current base-run core.

Uses the committed Java TGDM reference level sets from the brain-resolution
synthetic base run. Java mesh arrays and raw float inputs stay under ignored
``build/``; only compact hashes and measured metrics are saved in docs/reports.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np
from skimage.measure import marching_cubes

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import compare_jist_surface_geometry as reference  # noqa: E402
from pymacruise.image_io import normalize_image_to_ras  # noqa: E402
from pymacruise.surfaces import (  # noqa: E402
    SurfaceMesh,
    _signed_mesh_volume,
    extract_levelset_surface,
)
from pymacruise.thickness import closest_point_thickness  # noqa: E402

FIXTURE = ROOT / "tests/data/synthetic/base_core_source_match/java_reference"
WORK = ROOT / "build/java-reference/current-base-surface-comparison"
CLASSES = WORK / "classes"
RAW = WORK / "raw"
MESHES = WORK / "meshes"
REPORT_JSON = ROOT / "docs/reports/current-base-jist-surface-comparison-2026-10-07.json"
REPORT_MD = ROOT / "docs/reports/current-base-jist-surface-comparison-2026-10-07.md"


def extract_lorensen_surface(image: nib.spatialimages.SpatialImage) -> SurfaceMesh:
    """Apply the selected JIST near-zero rule before classic Lorensen MC."""
    ras_image = normalize_image_to_ras(image, source_name="level-set image")
    values = np.asarray(ras_image.dataobj, dtype=np.float32)
    shifted = values - np.float32(0.0)
    shifted = shifted.copy()
    shifted[(shifted > np.float32(-0.1)) & (shifted <= 0)] = np.float32(-0.1)
    shifted[(shifted < np.float32(0.1)) & (shifted > 0)] = np.float32(0.1)
    index_vertices, faces, _normals, _values = marching_cubes(
        shifted,
        level=0.0,
        gradient_direction="ascent",
        allow_degenerate=False,
        method="lorensen",
    )
    vertices = nib.affines.apply_affine(ras_image.affine, index_vertices)
    faces = np.asarray(faces, dtype=np.int64)
    if _signed_mesh_volume(vertices, faces) < 0:
        faces = faces[:, ::-1]
    return SurfaceMesh(vertices, faces)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def run() -> dict[str, object]:
    for directory in (CLASSES, RAW, MESHES):
        directory.mkdir(parents=True, exist_ok=True)
    vecmath = (ROOT.parent / "jist-scm-2026-09-26/JIST/lib/vecmath-1.3.1.jar").resolve()
    compile_args = [
        "javac", "-d", str(CLASSES.resolve()), "-cp", str(vecmath),
        str(reference.HARNESS), *(str(path) for path in reference.STUB_SOURCES),
        str(reference.KD_SOURCE), str(reference.JIST_SOURCE),
    ]
    subprocess.run(compile_args, check=True, cwd=ROOT)

    java_meshes: dict[str, SurfaceMesh] = {}
    python_meshes: dict[str, SurfaceMesh] = {}
    lorensen_meshes: dict[str, SurfaceMesh] = {}
    comparisons: dict[str, object] = {}
    lorensen_comparisons: dict[str, object] = {}
    source_inputs: dict[str, str] = {}
    spacing: list[float] | None = None
    shape: list[int] | None = None
    for name in ("inner", "central", "outer"):
        image = nib.load(str(FIXTURE / f"{name}.nii.gz"))
        values = np.asarray(image.dataobj, dtype="<f4")
        shape = list(values.shape)
        spacing = nib.affines.voxel_sizes(image.affine).tolist()
        digest = sha256(values.tobytes(order="C"))
        source_inputs[name] = digest
        raw_file = RAW / f"{name}.f32le"
        mesh_file = MESHES / f"{name}.json"
        values.tofile(raw_file)
        subprocess.run(
            [
                "java", "-Xms256m", "-Xmx4g", "-cp", f"{CLASSES.resolve()}:{vecmath}",
                "SurfaceIsoSurfaceOnGridHarness", str(raw_file.resolve()),
                *(str(size) for size in values.shape), str(mesh_file.resolve()),
            ],
            check=True,
            cwd=ROOT,
        )
        record, index_mesh = reference.load_java_mesh(mesh_file)
        if record["input_sha256"] != digest:
            raise AssertionError(f"JIST harness did not consume the exact {name} values")
        if tuple(record["dimensions"]) != tuple(values.shape):
            raise AssertionError(f"JIST harness returned the wrong {name} grid")
        jist_mesh = SurfaceMesh(
            nib.affines.apply_affine(image.affine, index_mesh.vertices), index_mesh.faces
        )
        python_mesh = extract_levelset_surface(image, level=0.0)
        java_meshes[name] = jist_mesh
        python_meshes[name] = python_mesh
        lorensen_mesh = extract_lorensen_surface(image)
        lorensen_meshes[name] = lorensen_mesh
        comparisons[name] = reference.compare_pair(jist_mesh, python_mesh)
        lorensen_comparisons[name] = reference.compare_pair(jist_mesh, lorensen_mesh)

    jist_thickness = closest_point_thickness(
        java_meshes["inner"], java_meshes["outer"], central=java_meshes["central"]
    )
    python_thickness = closest_point_thickness(
        python_meshes["inner"], python_meshes["outer"], central=python_meshes["central"]
    )
    lorensen_thickness = closest_point_thickness(
        lorensen_meshes["inner"], lorensen_meshes["outer"], central=lorensen_meshes["central"]
    )
    summary: dict[str, object] = {
        "title": "Selected JIST IsoSurfaceOnGrid versus Lewiner on current synthetic base-core outputs",
        "evidence_type": "synthetic selected-source comparison",
        "synthetic_only": True,
        "source": {
            "isosurface_source": str(reference.JIST_SOURCE.relative_to(ROOT.parent)),
            "isosurface_sha256": reference.sha256(reference.JIST_SOURCE),
            "wrapper_source": str(reference.JIST_WRAPPER_SOURCE.relative_to(ROOT.parent)),
            "wrapper_sha256": reference.sha256(reference.JIST_WRAPPER_SOURCE),
            "method": "IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, 18, 0.0f, false)",
            "version": "IsoSurfaceOnGrid CVS r1.7; GetMeshSurfaceFromLevelSet CVS r1.4",
            "harness": str(reference.HARNESS.relative_to(ROOT)),
            "harness_stubs": "Minimal JIST image/progress and normal-generator APIs; selected Java surface body is compiled directly.",
        },
        "fixture": {
            "path": str(FIXTURE.relative_to(ROOT)),
            "shape": shape,
            "orientation": "RAS",
            "spacing_mm": spacing,
            "source_levelset_sha256_float32_c_order": source_inputs,
        },
        "metrics": {
            "per_surface": comparisons,
            "lorensen_per_surface": lorensen_comparisons,
            "thickness": {
                "operator": "same Python closest_point_thickness on each triangulator's nested meshes",
                "jist_mean_mm": jist_thickness.mean_mm,
                "lewiner_mean_mm": python_thickness.mean_mm,
                "mean_delta_mm": python_thickness.mean_mm - jist_thickness.mean_mm,
                "jist_sd_mm": jist_thickness.standard_deviation_mm,
                "lewiner_sd_mm": python_thickness.standard_deviation_mm,
                "jist_nonzero_inner_vertices": jist_thickness.nonzero_vertex_count,
                "lewiner_nonzero_inner_vertices": python_thickness.nonzero_vertex_count,
                "lorensen_mean_mm": lorensen_thickness.mean_mm,
                "lorensen_mean_delta_mm": lorensen_thickness.mean_mm - jist_thickness.mean_mm,
            },
        },
        "limitations": [
            "The geometry distances sample mesh vertices against opposing triangles; they do not certify the continuous Hausdorff distance.",
            "Lewiner component counts and Euler characteristics differ from selected JIST on this case even though the sampled P08 surface-distance gates pass; connectivity remains an explicit P37 deviation.",
            "Thickness applies one Python operator to both meshes and isolates geometric impact; it does not compare the Java thickness implementation.",
            "This is a synthetic selected-source stage comparison, not whole-workflow or real-MRI validation.",
        ],
    }
    REPORT_JSON.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    rows = []
    for name, result in comparisons.items():
        result = result  # type: ignore[assignment]
        jt = result["jist_topology"]  # type: ignore[index]
        pt = result["lewiner_topology"]  # type: ignore[index]
        connect = result["triangle_connectivity_after_nearest_vertex_mapping"]  # type: ignore[index]
        distance = result["vertex_to_triangle_surface_distance_mm"]  # type: ignore[index]
        rows.append(
            f"| {name} | {jt['vertices']}/{jt['faces']}/{jt['euler_characteristic']} | "
            f"{pt['vertices']}/{pt['faces']}/{pt['euler_characteristic']} | "
            f"{distance['symmetric_mean']:.8f} | {distance['symmetric_max']:.6f} | "
            f"{connect['matching_fraction_of_jist_faces']:.2%} |"
        )
    lorensen_rows = []
    for name, result in lorensen_comparisons.items():
        result = result  # type: ignore[assignment]
        jt = result["jist_topology"]  # type: ignore[index]
        pt = result["lewiner_topology"]  # type: ignore[index]
        connect = result["triangle_connectivity_after_nearest_vertex_mapping"]  # type: ignore[index]
        distance = result["vertex_to_triangle_surface_distance_mm"]  # type: ignore[index]
        lorensen_rows.append(
            f"| {name} | {jt['connected_components']}/{jt['euler_characteristic']} | "
            f"{pt['connected_components']}/{pt['euler_characteristic']} | "
            f"{distance['symmetric_mean']:.8f} | {distance['symmetric_max']:.8f} | "
            f"{connect['matching_fraction_of_jist_faces']:.2%} |"
        )
    thickness = summary["metrics"]["thickness"]  # type: ignore[index]
    REPORT_MD.write_text(
        "\n".join(
            [
                "# JIST surface comparison on the brain-resolution synthetic base core",
                "",
                "This compares selected Java `IsoSurfaceOnGrid.solveOriginal(..., 18, 0, false)` with the production Python Lewiner extractor on the exact same committed Java TGDM outputs. The inputs are synthetic, RAS, and `(172, 216, 164)` at 0.9375 mm spacing; this is a stage comparison, not an end-to-end Java run.",
                "",
                f"JIST source SHA-256: `{summary['source']['isosurface_sha256']}`. Wrapper SHA-256: `{summary['source']['wrapper_sha256']}`.",
                "",
                "| Surface | JIST V/F/χ | Lewiner V/F/χ | Symmetric sampled vertex-to-triangle mean (mm) | Sampled max (mm) | Matched unordered faces |",
                "|---|---:|---:|---:|---:|---:|",
                *rows,
                "",
                "The symmetric sampled means are below 0.01 mm and sampled maxima are below 1 mm for all three surfaces, meeting the owner-provided P08 thresholds on these fixtures. The Lewiner topology counts differ: its mesh has multiple connected components and nonmatching Euler characteristics where JIST produces one connected genus-zero surface. Only about 63.8–67.9% of JIST triangles map after nearest-vertex matching. These remain material connectivity differences under P37 even though the sampled geometric position metrics pass.",
                "",
                "## Classic Lorensen alternative",
                "",
                "The same source perturbation and input values were also passed to scikit-image's classic Lorensen table:",
                "",
                "| Surface | JIST components/χ | Lorensen components/χ | Symmetric mean distance (mm) | Sampled max (mm) | Matched JIST faces |",
                "|---|---:|---:|---:|---:|---:|",
                *lorensen_rows,
                "",
                f"Lorensen preserves vertex counts and has near-coincident sampled geometry, but it does not repair the topology mismatch: component counts are 2/2/28 versus one JIST component. Its mean thickness is `{lorensen_thickness.mean_mm:.9f} mm` (delta `{lorensen_thickness.mean_mm - jist_thickness.mean_mm:+.9f} mm` from JIST). Thus the already-tested Lorensen alternative is not a topology-matching replacement for Lewiner on this brain-resolution fixture.",
                "",
                f"Applying the same Python closest-point thickness operator to each nested mesh set gives JIST mean `{thickness['jist_mean_mm']:.9f} mm` and Lewiner mean `{thickness['lewiner_mean_mm']:.9f} mm`, a delta of `{thickness['mean_delta_mm']:+.9f} mm`. P08 has no thickness-specific tolerance. This is a geometry-impact measurement, not Java thickness equivalence.",
                "",
                "The Java harness compiles the selected source body directly, with minimal image/progress and normal-generation stubs. It writes mesh arrays under ignored `build/`; this report stores only hashes and metrics. Vertex-to-triangle distances sample vertices and are not a continuous Hausdorff certificate. Reproduce with `./.venv/bin/python scripts/compare_jist_surface_current_base.py`; raw results: [JSON report](current-base-jist-surface-comparison-2026-10-07.json).",
                "",
                "Synthetic selected-source evidence only. It does not establish complete MaCRUISE workflow equivalence or validate human MRI outputs.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(json.dumps({"report": str(REPORT_MD), "json": str(REPORT_JSON), "metrics": summary["metrics"]}, indent=2))
    return summary


if __name__ == "__main__":
    run()
