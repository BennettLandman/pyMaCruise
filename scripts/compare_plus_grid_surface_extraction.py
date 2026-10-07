#!/usr/bin/env python3
"""Compare JIST and Lewiner on the saved 36-grid and 64-grid plus level sets.

Run the two natural-convergence synthetic workflows first, then:

    PYTHONPATH=src .venv/bin/python scripts/compare_plus_grid_surface_extraction.py

The selected Java IsoSurfaceOnGrid source is compiled into ignored build/
storage. The JSON result is also written under build/ so the source tree stays
small and no generated mesh reference is mistaken for committed source data.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.surfaces import SurfaceMesh, extract_levelset_surface
from pymacruise.validation import mesh_surface_distance_metrics, mesh_topology_metrics


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = ROOT.parent / "toads-cruise-scm-2026-09-26" / "toads-cruise" / "src"
JIST_SOURCE = SOURCE_ROOT / "edu/jhu/ece/iacl/algorithms/graphics/isosurf/IsoSurfaceOnGrid.java"
KD_SOURCE = SOURCE_ROOT / "edu/jhu/ece/iacl/algorithms/graphics/locator/kdtree/KdPoint3.java"
STUBS = ROOT / "tests/java_reference/surface_comparison_stubs"
HARNESS = ROOT / "tests/java_reference/SurfaceIsoSurfaceOnGridHarness.java"
JIST_LIB = ROOT.parent / "jist-scm-2026-09-26/JIST/lib"
VECMATH = JIST_LIB / "vecmath-1.3.1.jar"
CLASSES = ROOT / "build/java-reference/plus-grid-surface-classes"
OUTPUT = ROOT / "build/java-reference/plus-grid-surface-comparison"
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
RUN_MANIFESTS = {
    "36x45x34": ROOT / "build/plus-natural-convergence-complete/plus-surface/run.json",
    "64x80x64": ROOT / "build/plus-scaled-grid-current/natural-run-64/plus-surface/run.json",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _compile_reference() -> None:
    required = (JIST_SOURCE, KD_SOURCE, VECMATH, HARNESS, *STUB_SOURCES)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing selected JIST harness inputs: " + ", ".join(missing))
    CLASSES.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "javac", "-d", str(CLASSES), "-cp", str(VECMATH),
            str(HARNESS), *(str(path) for path in STUB_SOURCES),
            str(KD_SOURCE), str(JIST_SOURCE),
        ],
        check=True,
        cwd=ROOT,
    )


def compare() -> dict[str, object]:
    _compile_reference()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    cases: dict[str, object] = {}
    for grid_name, manifest_path in RUN_MANIFESTS.items():
        if not manifest_path.is_file():
            raise FileNotFoundError(
                f"Missing completed plus surface manifest: {manifest_path}; run the synthetic workflow first"
            )
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("status") != "complete":
            raise ValueError(f"Plus surface run is not complete: {manifest_path}")
        stage_outputs = manifest["stages"]["surfaces_thickness"]["outputs"]
        grid_result: dict[str, object] = {}
        for surface_name in ("inner", "central", "outer"):
            image_path = Path(stage_outputs[f"{surface_name}_level_set"]["path"])
            if _sha256(image_path) != stage_outputs[f"{surface_name}_level_set"]["sha256"]:
                raise ValueError(f"Level-set manifest hash mismatch: {image_path}")
            image = nib.load(str(image_path))
            if len(image.shape) != 3 or nib.aff2axcodes(image.affine) != ("R", "A", "S"):
                raise ValueError(f"Expected a 3D RAS level set: {image_path}")
            values = np.asarray(image.dataobj, dtype="<f4")
            if not np.isfinite(values).all():
                raise ValueError(f"Nonfinite level set: {image_path}")
            raw_path = OUTPUT / f"{grid_name}-{surface_name}.raw.f32le"
            java_path = OUTPUT / f"{grid_name}-{surface_name}.json"
            values.tofile(raw_path)
            subprocess.run(
                [
                    "java", "-Xmx512m", "-cp", f"{CLASSES}:{VECMATH}",
                    "SurfaceIsoSurfaceOnGridHarness", str(raw_path),
                    *(str(size) for size in values.shape), str(java_path),
                ],
                check=True,
                cwd=ROOT,
            )
            java_record = json.loads(java_path.read_text(encoding="utf-8"))
            input_sha256 = hashlib.sha256(values.tobytes(order="C")).hexdigest()
            if java_record.get("input_sha256") != input_sha256:
                raise ValueError(f"Java harness consumed different float32 data: {image_path}")
            if (
                java_record.get("dimensions") != list(values.shape)
                or java_record.get("connectivity") != 18
                or java_record.get("iso_level") != 0
                or java_record.get("inclusive") is not False
            ):
                raise ValueError(f"Java harness settings or dimensions differ: {java_path}")

            jist_mesh = SurfaceMesh(
                np.asarray(java_record["vertices"], dtype=np.float64),
                np.asarray(java_record["faces"], dtype=np.int64),
            )
            # The harness exposes source voxel coordinates. Apply the exact
            # NIfTI affine so both results are compared in required RAS mm.
            jist_mesh = SurfaceMesh(
                jist_mesh.vertices @ image.affine[:3, :3].T + image.affine[:3, 3],
                jist_mesh.faces,
            )
            python_mesh = extract_levelset_surface(image, level=0.0)
            grid_result[surface_name] = {
                "level_set_path": str(image_path.relative_to(ROOT)),
                "level_set_sha256": _sha256(image_path),
                "float32_data_sha256": input_sha256,
                "shape": list(values.shape),
                "jist": mesh_topology_metrics(jist_mesh.vertices, jist_mesh.faces),
                "lewiner": mesh_topology_metrics(python_mesh.vertices, python_mesh.faces),
                "symmetric_ras_mm": mesh_surface_distance_metrics(
                    jist_mesh.vertices,
                    jist_mesh.faces,
                    python_mesh.vertices,
                    python_mesh.faces,
                ),
            }
            print(grid_name, surface_name, grid_result[surface_name]["symmetric_ras_mm"])
        cases[grid_name] = grid_result

    return {
        "method": "Selected IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, connectivity=18, iso=0, inclusive=false) versus scikit-image Lewiner on exact same float32 arrays",
        "coordinate_comparison": "Java harness voxel-index vertices are transformed by each level-set NIfTI affine; metrics are symmetric vertex-to-continuous-triangle distances in RAS mm.",
        "source": {
            "jist_source_sha256": _sha256(JIST_SOURCE),
            "jist_source_revision": "CVS r1.7",
            "harness_source": str(HARNESS.relative_to(ROOT)),
        },
        "fixture_runs": {
            name: {
                "manifest": str(path.relative_to(ROOT)),
                "manifest_sha256": _sha256(path),
            }
            for name, path in RUN_MANIFESTS.items()
        },
        "cases": cases,
        "limitations": [
            "The harness executes the selected JIST surface algorithm with minimal image/progress stubs; it does not compare the complete JIST runtime or VTK writer.",
            "This measures the surface-extraction substitution on two synthetic workflow level sets; it does not establish whole-workflow Java/MATLAB equivalence.",
            "The vertex-to-triangle metric samples mesh vertices and is not a continuous triangle-to-triangle Hausdorff certificate.",
        ],
    }


def main() -> int:
    result = compare()
    result_path = OUTPUT / "comparison.json"
    result_path.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {result_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
