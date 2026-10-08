#!/usr/bin/env python3
"""Compare selected JIST meshes with Python meshes from a saved workflow.

This operates on the exact saved inner/central/outer scalar fields in
``--core-dir``. It compiles the selected Java IsoSurfaceOnGrid source harness,
then compares each Java mesh with the production scikit-image Lewiner mesh.
Temporary raw arrays and mesh outputs stay under ``build/``.
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
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import compare_jist_surface_geometry as reference  # noqa: E402
from pymacruise.jist_intersector import JistSurfaceTree  # noqa: E402
from pymacruise.surfaces import extract_levelset_surface  # noqa: E402


def compare_meshes(reference_mesh, candidate_mesh, *, max_samples: int = 4096) -> dict[str, object]:
    """Compare topology exactly and nearest-surface distances on fixed samples."""
    reference_count = len(reference_mesh.vertices)
    candidate_count = len(candidate_mesh.vertices)
    reference_samples = np.linspace(
        0, reference_count - 1, min(reference_count, max_samples), dtype=np.int64
    )
    candidate_samples = np.linspace(
        0, candidate_count - 1, min(candidate_count, max_samples), dtype=np.int64
    )
    reference_to_candidate = JistSurfaceTree(
        candidate_mesh.vertices, candidate_mesh.faces
    ).query_many(reference_mesh.vertices[reference_samples])[1]
    candidate_to_reference = JistSurfaceTree(
        reference_mesh.vertices, reference_mesh.faces
    ).query_many(candidate_mesh.vertices[candidate_samples])[1]

    nearest_indices = cKDTree(candidate_mesh.vertices).query(reference_mesh.vertices)[1]
    mapped_reference_faces = {
        tuple(sorted(int(nearest_indices[vertex]) for vertex in face))
        for face in reference_mesh.faces
    }
    candidate_faces = {
        tuple(sorted(int(vertex) for vertex in face)) for face in candidate_mesh.faces
    }
    matching_faces = len(mapped_reference_faces & candidate_faces)
    all_distances = np.concatenate((reference_to_candidate, candidate_to_reference))
    return {
        "jist_topology": reference.mesh_topology(reference_mesh),
        "python_topology": reference.mesh_topology(candidate_mesh),
        "triangle_connectivity_after_nearest_vertex_mapping": {
            "matching_unordered_faces": matching_faces,
            "jist_face_count": int(len(reference_mesh.faces)),
            "python_face_count": int(len(candidate_mesh.faces)),
            "matching_fraction_of_jist_faces": float(
                matching_faces / len(reference_mesh.faces)
            ),
        },
        "sampled_symmetric_vertex_to_triangle_distance_mm": {
            "sample_count_per_direction": int(min(reference_count, max_samples)),
            "symmetric_mean": float(np.mean(all_distances)),
            "symmetric_max": float(np.max(all_distances)),
        },
    }


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--core-dir", type=Path, required=True,
        help="saved workflow CRUISE core stage containing inner/central/outer NIfTIs",
    )
    parser.add_argument(
        "--work-dir", type=Path, required=True,
        help="ignored build directory for raw inputs and Java mesh arrays",
    )
    args = parser.parse_args()
    core_dir = args.core_dir.resolve()
    work_dir = args.work_dir.resolve()
    classes = work_dir / "classes"
    raw_dir = work_dir / "raw"
    meshes_dir = work_dir / "meshes"
    for path in (classes, raw_dir, meshes_dir):
        path.mkdir(parents=True, exist_ok=True)

    vecmath = ROOT.parent / "jist-scm-2026-09-26/JIST/lib/vecmath-1.3.1.jar"
    required = (reference.JIST_SOURCE, reference.KD_SOURCE, reference.HARNESS, vecmath)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing selected JIST source/runtime: " + ", ".join(missing))
    compile_command = [
        "javac", "-d", str(classes), "-cp", str(vecmath),
        str(reference.HARNESS), *(str(path) for path in reference.STUB_SOURCES),
        str(reference.KD_SOURCE), str(reference.JIST_SOURCE),
    ]
    subprocess.run(compile_command, check=True, cwd=ROOT)

    metrics: dict[str, object] = {}
    input_hashes: dict[str, str] = {}
    shape: list[int] | None = None
    spacing: list[float] | None = None
    for name in ("inner", "central", "outer"):
        image_path = core_dir / f"{name}.nii.gz"
        image = nib.load(str(image_path))
        if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise ValueError(f"{image_path} is not RAS-oriented")
        values = np.asarray(image.dataobj, dtype="<f4")
        if values.ndim != 3 or not np.isfinite(values).all():
            raise ValueError(f"{image_path} must contain finite 3-D values")
        shape = list(values.shape)
        spacing = nib.affines.voxel_sizes(image.affine).tolist()
        raw_path = raw_dir / f"{name}.f32le"
        mesh_path = meshes_dir / f"{name}.json"
        raw_bytes = values.tobytes(order="C")
        input_hashes[name] = sha256_bytes(raw_bytes)
        raw_path.write_bytes(raw_bytes)
        subprocess.run(
            [
                "java", "-Xms256m", "-Xmx4g", "-cp", f"{classes}:{vecmath}",
                "SurfaceIsoSurfaceOnGridHarness", str(raw_path),
                *(str(size) for size in values.shape), str(mesh_path),
            ],
            check=True,
            cwd=ROOT,
        )
        record, index_mesh = reference.load_java_mesh(mesh_path)
        if record["input_sha256"] != input_hashes[name]:
            raise AssertionError(f"JIST did not consume the exact {name} field")
        if tuple(record["dimensions"]) != tuple(values.shape):
            raise AssertionError(f"JIST returned the wrong {name} grid")
        jist_mesh = reference.SurfaceMesh(
            nib.affines.apply_affine(image.affine, index_mesh.vertices), index_mesh.faces
        )
        python_mesh = extract_levelset_surface(image, level=0.0)
        metrics[name] = compare_meshes(jist_mesh, python_mesh)

    result = {
        "evidence_type": "synthetic selected-source comparison",
        "synthetic_only": True,
        "core_dir": str(core_dir),
        "shape": shape,
        "spacing_mm": spacing,
        "orientation": "RAS",
        "source_scalar_sha256_float32_c_order": input_hashes,
        "jist_source_sha256": reference.sha256(reference.JIST_SOURCE),
        "harness_sha256": reference.sha256(reference.HARNESS),
        "source_method": "IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, 18, 0.0f, false)",
        "comparison": metrics,
        "limits": [
            "Vertex-to-triangle distances sample mesh vertices; they do not certify continuous Hausdorff distance.",
            "This is a held-fixed surface-extraction comparison, not whole-workflow equivalence.",
        ],
    }
    report_path = work_dir / "comparison.json"
    report_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"report": str(report_path), "results": result}, indent=2))


if __name__ == "__main__":
    main()
