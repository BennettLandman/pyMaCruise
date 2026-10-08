#!/usr/bin/env python3
"""Measure scikit-image Lorensen marching cubes against saved JIST meshes.

The Java meshes were produced from the same folded RAS level sets by
``IsoSurfaceOnGrid.solveOriginal``. This diagnostic asks whether scikit-image's
classic case table removes the face-connectivity difference seen with its
Lewiner implementation; it does not alter the production surface extractor.
"""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np
from skimage.measure import marching_cubes

from compare_jist_surface_geometry import (
    FIXTURES,
    JAVA_MESH_DIR,
    compare_pair,
    load_java_mesh,
)
from pymacruise.surfaces import SurfaceMesh
from pymacruise.thickness import closest_point_thickness


def extract_lorensen(image: nib.spatialimages.SpatialImage) -> SurfaceMesh:
    values = np.asarray(image.dataobj, dtype=np.float32)
    shifted = values - np.float32(0.0)
    shifted = shifted.copy()
    shifted[(shifted > np.float32(-0.1)) & (shifted <= 0)] = np.float32(-0.1)
    shifted[(shifted < np.float32(0.1)) & (shifted > 0)] = np.float32(0.1)
    vertices, faces, _, _ = marching_cubes(
        shifted,
        level=0.0,
        gradient_direction="ascent",
        allow_degenerate=False,
        method="lorensen",
    )
    vertices = vertices @ image.affine[:3, :3].T + image.affine[:3, 3]
    return SurfaceMesh(vertices, np.asarray(faces, dtype=np.int64))


def main() -> None:
    references: dict[str, SurfaceMesh] = {}
    candidates: dict[str, SurfaceMesh] = {}
    surfaces: dict[str, dict[str, object]] = {}
    for name in ("inner", "central", "outer"):
        image = nib.load(str(FIXTURES / f"{name}.nii.gz"))
        record, reference = load_java_mesh(JAVA_MESH_DIR / f"{name}.json")
        values = np.asarray(image.dataobj, dtype="<f4")
        input_hash = hashlib.sha256(values.tobytes(order="C")).hexdigest()
        if record["input_sha256"] != input_hash:
            raise AssertionError(f"Java harness did not consume the exact {name} float32 input")
        if tuple(record["dimensions"]) != tuple(values.shape):
            raise AssertionError(f"Java harness dimension mismatch for {name}")
        candidate = extract_lorensen(image)
        references[name] = reference
        candidates[name] = candidate
        surfaces[name] = {
            "input_sha256": input_hash,
            **compare_pair(reference, candidate),
        }

    reference_thickness = closest_point_thickness(
        references["inner"], references["outer"], central=references["central"]
    )
    candidate_thickness = closest_point_thickness(
        candidates["inner"], candidates["outer"], central=candidates["central"]
    )
    output = {
        "method": "scikit-image marching_cubes(method='lorensen')",
        "comparison": "selected JIST IsoSurfaceOnGrid.solveOriginal vs Lorensen on identical folded RAS level sets",
        "fixture": "tests/data/synthetic/jist_surface_comparison; 41x41x41 folded sphere",
        "isosurface_java_meshes": "tests/data/synthetic/jist_surface_comparison/java_meshes",
        "per_surface": surfaces,
        "thickness": {
            "operator": "same Python closest_point_thickness applied to both mesh sets",
            "jist_mean_mm": float(reference_thickness.mean_mm),
            "lorensen_mean_mm": float(candidate_thickness.mean_mm),
            "mean_delta_mm": float(candidate_thickness.mean_mm - reference_thickness.mean_mm),
        },
        "limitations": [
            "Synthetic source-comparison diagnostic only; no human MRI data.",
            "Sampled vertex-to-triangle distances are not a continuous Hausdorff certificate.",
            "The Lorensen alternative is not an exact JIST connectivity replacement and does not close P37.",
        ],
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
