#!/usr/bin/env python3
"""Compare the JIST reference with the production Python port on one checkpoint."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import compare_jist_surface_geometry as reference  # noqa: E402
from pymacruise.surfaces import SurfaceMesh, extract_levelset_surface  # noqa: E402
from pymacruise.thickness import closest_point_thickness  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="terminal or checkpoint surface manifest")
    parser.add_argument("--output", type=Path, required=True, help="JSON report path")
    parser.add_argument("--markdown", type=Path, required=True, help="Markdown report path")
    args = parser.parse_args()
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    source_paths = {
        "inner": Path(manifest["final_outputs"]["plus_inner_level_set"]["path"]),
        "central": Path(manifest["final_outputs"]["plus_central_level_set"]["path"]),
        "outer": Path(manifest["final_outputs"]["plus_outer_level_set"]["path"]),
    }
    fixture_dir = ROOT / "build/java-reference/terminal-plus-levelsets"
    fixture_dir.mkdir(parents=True, exist_ok=True)
    for name, source in source_paths.items():
        destination = fixture_dir / f"{name}.nii.gz"
        destination.unlink(missing_ok=True)
        destination.symlink_to(source.resolve())

    work = ROOT / "build/java-reference/terminal-plus-surface-comparison"
    reference.FIXTURES = fixture_dir
    reference.CLASSES = work / "classes"
    reference.JAVA_MESH_DIR = work / "meshes"
    reference.RAW_DIR = work / "raw"
    java_outputs = reference.run_java()

    meshes: dict[str, dict[str, SurfaceMesh]] = {}
    details: dict[str, dict[str, object]] = {}
    input_hashes: dict[str, str] = {}
    for name, path in source_paths.items():
        image = nib.load(str(path))
        values = np.asarray(image.dataobj, dtype="<f4")
        digest = hashlib.sha256(values.tobytes(order="C")).hexdigest()
        record, index_mesh = reference.load_java_mesh(java_outputs[name])
        if record["input_sha256"] != digest:
            raise AssertionError(f"JIST did not consume the exact {name} level set")
        if tuple(record["dimensions"]) != tuple(values.shape):
            raise AssertionError(f"JIST dimensions mismatch for {name}")
        jist = SurfaceMesh(nib.affines.apply_affine(image.affine, index_mesh.vertices), index_mesh.faces)
        python_mesh = extract_levelset_surface(image, level=0.0)
        meshes[name] = {"jist": jist, "python_source_matched": python_mesh}
        comparison = reference.compare_pair(jist, python_mesh)
        distance = comparison["vertex_to_triangle_surface_distance_mm"]
        connectivity = comparison["triangle_connectivity_after_nearest_vertex_mapping"]
        details[name] = {
            "jist_topology": comparison["jist_topology"],
            "python_source_matched_topology": comparison["lewiner_topology"],
            "triangle_connectivity_after_nearest_vertex_mapping": {
                "matching_unordered_faces": connectivity["matching_unordered_faces"],
                "jist_face_count": connectivity["jist_face_count"],
                "python_face_count": connectivity["lewiner_face_count"],
                "matching_fraction_of_jist_faces": connectivity[
                    "matching_fraction_of_jist_faces"
                ],
                "vertex_mapping_is_bijective": connectivity["vertex_mapping_is_bijective"],
            },
            "vertex_to_triangle_surface_distance_mm": {
                "symmetric_mean": distance["symmetric_mean"],
                "symmetric_max": distance["symmetric_max"],
                "jist_vertices_to_python_mean": distance["jist_vertices_to_lewiner_mean"],
                "python_vertices_to_jist_mean": distance["lewiner_vertices_to_jist_mean"],
            },
        }
        input_hashes[name] = digest

    thickness = {}
    for backend in ("jist", "python_source_matched"):
        value = closest_point_thickness(
            meshes["inner"][backend], meshes["outer"][backend],
            central=meshes["central"][backend],
        )
        thickness[backend] = {
            "mean_mm": value.mean_mm,
            "standard_deviation_mm": value.standard_deviation_mm,
            "nonzero_inner_vertices": value.nonzero_vertex_count,
        }
    report = {
        "evidence_category": "synthetic selected-source surface-extraction comparison",
        "workflow_manifest": str(manifest_path),
        "workflow_terminal": manifest.get("terminal"),
        "completed_iteration": manifest.get("completed_iteration"),
        "source_method": "IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, 18, 0.0f, false)",
        "python_method": "production source-derived 18/6 templates generated from all 256 selected Java corner-sign cases",
        "source_sha256": reference.sha256(reference.JIST_SOURCE),
        "wrapper_sha256": reference.sha256(reference.JIST_WRAPPER_SOURCE),
        "same_levelset_float32_sha256": input_hashes,
        "per_surface": details,
        "thickness_same_python_operator_mm": thickness,
        "thickness_mean_delta_python_minus_jist_mm": thickness["python_source_matched"]["mean_mm"] - thickness["jist"]["mean_mm"],
        "limitations": [
            "This isolates triangulation on the exact same level-set arrays; it does not compare level-set evolution or the full Java workflow.",
            "Reported symmetric distances sample vertices against triangles and do not certify continuous Hausdorff distance.",
            "The Python thickness function is applied to both mesh sets; this isolates mesh geometry and is not Java thickness equivalence.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    rows = []
    for name, comparison in details.items():
        jt, lw = comparison["jist_topology"], comparison["python_source_matched_topology"]
        dist = comparison["vertex_to_triangle_surface_distance_mm"]
        rows.append(
            f"| {name} | {jt['vertices']}/{jt['faces']}/{jt['connected_components']}/{jt['euler_characteristic']} | "
            f"{lw['vertices']}/{lw['faces']}/{lw['connected_components']}/{lw['euler_characteristic']} | "
            f"{dist['symmetric_mean']:.8f} | {dist['symmetric_max']:.6f} |"
        )
    args.markdown.write_text(
        "\n".join([
            "# Terminal synthetic PLUS surface extraction: Java JIST and Python source port",
            "",
            "The source-selected `IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, 18, 0, false)` and the production Python source-derived 18/6 case-template extractor were given the exact same terminal checkpoint level-set arrays. This is a synthetic stage comparison of the current source port; it does not establish upstream or whole-workflow equivalence.",
            "",
            f"The registered LesionTOADS loop stopped by `{manifest.get('stop_reason')}` at iteration {manifest.get('completed_iteration')}. The input level sets were checked by SHA-256 before comparing the generated meshes.",
            "",
            "| Surface | Java JIST V/F/components/χ | Python source-port V/F/components/χ | Symmetric sampled mean (mm) | Sampled max (mm) |",
            "|---|---:|---:|---:|---:|",
            *rows,
            "",
            f"Using the same Python closest-point thickness operator on both mesh sets, mean thickness is {thickness['jist']['mean_mm']:.9f} mm (JIST) and {thickness['python_source_matched']['mean_mm']:.9f} mm (Python), delta {report['thickness_mean_delta_python_minus_jist_mm']:+.9f} mm. This is not a Java thickness comparison.",
            "",
            "Distance summaries sample mesh vertices against triangles; their maxima are not continuous Hausdorff certificates. See the adjacent JSON for hashes, topology, and method provenance.",
            "",
        ]),
        encoding="utf-8",
    )
    print(json.dumps({"json": str(args.output), "markdown": str(args.markdown), "per_surface": details, "thickness": thickness}, indent=2))


if __name__ == "__main__":
    main()
