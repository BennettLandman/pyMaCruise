#!/usr/bin/env python3
"""Measure the production Lewiner-mesh thickness versus selected Java meshes.

This reuses the full-grid Java wrapper output from
``compare_jist_closest_point_thickness_brain_resolution.py`` and computes
thickness only on the current production surface extractor's meshes.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import nibabel as nib
import numpy as np

from pymacruise.surfaces import extract_levelset_surface
from pymacruise.thickness import closest_point_thickness
import compare_jist_closest_point_thickness as selected


ROOT = Path(__file__).resolve().parents[1]
LEVELSETS = ROOT / "tests/data/synthetic/base_core_source_match/java_reference"
JAVA_MESHES = ROOT / "build/java-reference/current-base-surface-comparison/meshes"
JAVA_OUTPUT = selected.WORK / "base-result.raw"
REPORT = ROOT / "docs/reports/current-base-thickness-substitution-2026-10-07.json"
REPORT_MD = ROOT / "docs/reports/current-base-thickness-substitution-2026-10-07.md"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    if not JAVA_OUTPUT.is_file():
        raise FileNotFoundError(
            "selected Java thickness output is missing; first run "
            "scripts/compare_jist_closest_point_thickness_brain_resolution.py"
        )
    raw = JAVA_OUTPUT.read_bytes()
    inner_count = int(np.frombuffer(raw, dtype="<i4", count=1, offset=0)[0])
    java_inner = np.frombuffer(raw, dtype="<f8", count=inner_count * 3, offset=4).reshape(inner_count, 3)
    java_vertex_thickness = java_inner[:, 2]
    java_nonzero = java_vertex_thickness[java_vertex_thickness != 0]
    java_mean = float(java_nonzero.mean())
    java_std = float(java_nonzero.std(ddof=1)) if java_nonzero.size > 1 else 0.0

    images = {name: nib.load(str(LEVELSETS / f"{name}.nii.gz")) for name in ("inner", "central", "outer")}
    production_meshes = {name: extract_levelset_surface(image, level=0.0) for name, image in images.items()}
    started = time.perf_counter()
    result = closest_point_thickness(
        production_meshes["inner"], production_meshes["outer"], central=production_meshes["central"]
    )
    elapsed = time.perf_counter() - started
    peak_rss_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (
        1e9 if sys.platform == "darwin" else 1e6
    )
    record = {
        "title": "Brain-resolution thickness impact of selected surface-extraction substitution",
        "synthetic_only": True,
        "comparison_scope": "The selected Java ClosestPointThicknessForInner wrapper output on Java IsoSurfaceOnGrid meshes is compared by aggregate inner mean/stdev with the production Python thickness API on Lewiner meshes extracted from the same committed level sets. Vertex correspondence is not assumed because topology and face connectivity differ.",
        "fixture": {
            "levelset_path": str(LEVELSETS.relative_to(ROOT)),
            "levelset_sha256": {name: sha(LEVELSETS / f"{name}.nii.gz") for name in images},
            "java_mesh_path": str(JAVA_MESHES.relative_to(ROOT)),
            "java_mesh_sha256": {name: sha(JAVA_MESHES / f"{name}.json") for name in images},
            "selected_java_output_sha256": sha(JAVA_OUTPUT),
            "production_extractor": "skimage.measure.marching_cubes method=lewiner via extract_levelset_surface",
            "mesh_counts": {
                name: {"vertices": len(mesh.vertices), "faces": len(mesh.faces)}
                for name, mesh in production_meshes.items()
            },
        },
        "metrics": {
            "selected_java_mean_thickness_mm": java_mean,
            "production_python_mean_thickness_mm": result.mean_mm,
            "mean_difference_python_minus_java_mm": result.mean_mm - java_mean,
            "selected_java_sample_stdev_mm": java_std,
            "production_python_sample_stdev_mm": result.standard_deviation_mm,
            "stdev_difference_python_minus_java_mm": result.standard_deviation_mm - java_std,
            "selected_java_nonzero_inner_vertices": int(java_nonzero.size),
            "production_python_nonzero_inner_vertices": result.nonzero_vertex_count,
            "thickness_specific_acceptance_threshold_defined": False,
        },
        "runtime": {
            "production_python_api_elapsed_seconds": elapsed,
            "python_process_peak_rss_gb": peak_rss_gb,
        },
        "limitations": [
            "This is an aggregate measurement comparison; it does not establish per-vertex correspondence across different meshes.",
            "The user-defined P08 tolerances cover segmentation Dice and surface position, not thickness scalar differences.",
            "The Java wrapper result is source-equivalent on its held-fixed Java meshes; Lewiner face connectivity remains a documented substitution.",
            "Synthetic stage evidence only; whole-workflow and real-MRI equivalence remain open.",
        ],
    }
    REPORT.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [
        "# Brain-resolution thickness impact of surface extraction",
        "",
        "This compares the selected Java thickness-wrapper aggregate over Java IsoSurfaceOnGrid meshes with the production Python thickness API over Lewiner meshes extracted from the same committed level sets. Mesh topology differs, so this report compares aggregate values and does not match vertices one-to-one.",
        "",
        "| Metric | Selected Java / JIST meshes | Production Python / Lewiner meshes | Difference |",
        "|---|---:|---:|---:|",
        f"| Mean thickness (mm) | {java_mean:.9f} | {result.mean_mm:.9f} | {result.mean_mm - java_mean:+.9f} |",
        f"| Sample standard deviation (mm) | {java_std:.9f} | {result.standard_deviation_mm:.9f} | {result.standard_deviation_mm - java_std:+.9f} |",
        f"| Nonzero inner vertices | {java_nonzero.size:,} | {result.nonzero_vertex_count:,} | — |",
        "",
        f"The production Python thickness call took `{elapsed:.1f} s` and peaked at `{peak_rss_gb:.3f} GB` RSS. No thickness-specific acceptance threshold has been set; the report quantifies the current surface-extraction substitution rather than declaring it equivalent. The P08 surface-position checks and topology differences are documented in the [current-base surface report](current-base-jist-surface-comparison-2026-10-07.md). Reproduce after generating the selected Java comparison with `./.venv/bin/python scripts/compare_jist_closest_point_thickness_brain_resolution.py` by running `./.venv/bin/python scripts/compare_current_base_thickness_substitution.py`.",
        "",
    ]
    REPORT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"report": str(REPORT), "metrics": record["metrics"], "runtime": record["runtime"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
