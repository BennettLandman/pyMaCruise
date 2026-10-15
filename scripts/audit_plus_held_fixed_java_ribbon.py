#!/usr/bin/env python3
"""Audit tissue-to-ribbon diagnostics on held-fixed PLUS Java and Python outputs.

This intentionally uses Java TGDM outputs generated from Python's exact crop
inputs. It answers whether a poor anatomy diagnostic is also present in that
selected-source replay; it is not an independent full PLUS source workflow.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pymacruise.cortical_separation import GM_LABELS, WM_LABELS  # noqa: E402
from pymacruise.surfaces import extract_levelset_surface  # noqa: E402
from pymacruise.thickness import closest_point_thickness, levelset_subtraction_thickness  # noqa: E402


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_array(path: Path, shape: tuple[int, int, int]) -> np.ndarray:
    with gzip.open(path, "rb") as stream:
        raw = stream.read()
    expected_bytes = int(np.prod(shape)) * np.dtype(">f4").itemsize
    if len(raw) != expected_bytes:
        raise ValueError(f"{path} has {len(raw)} bytes; expected {expected_bytes}")
    return np.frombuffer(raw, dtype=">f4").astype(np.float32).reshape(shape)


def _topology(mesh) -> dict[str, int | bool]:
    edges = np.concatenate(
        (mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]])
    )
    edges.sort(axis=1)
    unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
    adjacency = coo_matrix(
        (np.ones(len(unique_edges), dtype=np.uint8),
         (unique_edges[:, 0], unique_edges[:, 1])),
        shape=(len(mesh.vertices), len(mesh.vertices)),
    )
    components, _ = connected_components(adjacency + adjacency.T, directed=False)
    return {
        "vertices": int(len(mesh.vertices)),
        "faces": int(len(mesh.faces)),
        "connected_components": int(components),
        "boundary_edges": int(np.count_nonzero(counts == 1)),
        "nonmanifold_edges": int(np.count_nonzero(counts > 2)),
        "euler_characteristic": int(len(mesh.vertices) - len(unique_edges) + len(mesh.faces)),
        "watertight": bool(np.all(counts == 2)),
    }


def _audit(arrays: dict[str, np.ndarray], labels: np.ndarray, affine: np.ndarray) -> dict:
    images = {name: nib.Nifti1Image(values, affine) for name, values in arrays.items()}
    masks = {name: values <= 0 for name, values in arrays.items()}
    wm = np.isin(labels, tuple(WM_LABELS))
    gm = np.isin(labels, tuple(GM_LABELS))
    n_wm, n_gm = int(wm.sum()), int(gm.sum())
    support = {
        "wm_inside_inner": float(np.count_nonzero(wm & masks["inner"]) / max(1, n_wm)),
        "gm_between_inner_and_outer": float(
            np.count_nonzero(gm & ~masks["inner"] & masks["outer"]) / max(1, n_gm)
        ),
    }
    nesting = {
        "inner_outside_central_voxels": int(np.count_nonzero(masks["inner"] & ~masks["central"])),
        "central_outside_outer_voxels": int(np.count_nonzero(masks["central"] & ~masks["outer"])),
    }
    meshes = {name: extract_levelset_surface(images[name]) for name in ("inner", "central", "outer")}
    thickness = closest_point_thickness(meshes["inner"], meshes["outer"], central=meshes["central"])
    volume = levelset_subtraction_thickness(images["inner"], images["outer"])
    values = thickness.inner_surface.point_data["mean_thickness_mm"]
    return {
        "label_counts": {"wm": n_wm, "gm": n_gm},
        "support_fractions": support,
        "nesting_voxel_counts": nesting,
        "surfaces": {name: _topology(mesh) for name, mesh in meshes.items()},
        "surface_thickness_mm": {
            "mean": float(thickness.mean_mm),
            "median": float(np.median(values)),
            "p05": float(np.percentile(values, 5)),
            "p95": float(np.percentile(values, 95)),
            "sample_count": int(thickness.nonzero_vertex_count),
        },
        "levelset_thickness_mm": {
            "mean": float(volume.mean_mm),
            "crossing_voxels": int(volume.crossing_voxel_count),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    reference = args.reference_dir.resolve()
    provenance_path = args.provenance.resolve()
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    shape = tuple(map(int, provenance["shape_xyz"]))
    affine = np.asarray(provenance["affine_ras"], dtype=np.float64)

    segmentation = provenance["crop_inputs"]["segmentation"]
    segmentation_path = Path(segmentation["path"]).resolve()
    if _sha256(segmentation_path) != segmentation["sha256"]:
        raise ValueError("crop MaSeg input does not match its pinned SHA-256")
    segmentation_image = nib.load(str(segmentation_path))
    if segmentation_image.shape != shape or not np.allclose(segmentation_image.affine, affine, rtol=0, atol=1e-5):
        raise ValueError("crop MaSeg geometry differs from the Java source grid")
    if nib.aff2axcodes(segmentation_image.affine) != ("R", "A", "S"):
        raise ValueError("crop MaSeg is not RAS")
    labels = np.asarray(segmentation_image.dataobj)
    if not np.equal(labels, np.floor(labels)).all():
        raise ValueError("crop MaSeg contains noninteger labels")

    source_arrays = {
        name: _source_array(reference / "java" / f"{name}.raw.gz", shape)
        for name in ("inner", "central", "outer")
    }
    python_arrays = {
        name: np.load(reference / f"python-{name}.npy", mmap_mode="r").astype(np.float32)
        for name in ("inner", "central", "outer")
    }
    for name, values in {**source_arrays, **python_arrays}.items():
        if values.shape != shape or not np.isfinite(values).all():
            raise ValueError(f"{name} field has invalid shape or nonfinite values")

    report = {
        "case": "case-100260 human PLUS held-fixed Java ribbon audit",
        "evidence_category": "real-human anatomy diagnostics on selected Java TGDM replay with Python workflow inputs held fixed",
        "scope_limit": "not an independent full PLUS Java/MATLAB workflow; Java TGDM starts from the same crop inputs produced by the Python run",
        "coordinate_system": "RAS",
        "shape_xyz": list(shape),
        "affine_ras": affine.tolist(),
        "provenance_sha256": _sha256(provenance_path),
        "crop_maseg_sha256": segmentation["sha256"],
        "source_tgdm": _audit(source_arrays, labels, affine),
        "python_tgdm": _audit(python_arrays, labels, affine),
        "field_comparison": {
            name: {
                "differing_voxels": int(np.count_nonzero(source_arrays[name] != python_arrays[name])),
                "sign_differing_voxels": int(np.count_nonzero(np.signbit(source_arrays[name]) != np.signbit(python_arrays[name]))),
            }
            for name in source_arrays
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
