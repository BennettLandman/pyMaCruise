#!/usr/bin/env python3
"""Audit a completed real-data base MaCRUISE continuation.

The audit records file hashes and image geometry, checks finite RAS NIfTIs,
measures cortical MaSeg support and nested level-set regions, and reports mesh
topology from the already-written VTK products. It does not re-extract surfaces
or claim source equivalence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import nibabel as nib
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from pymacruise.cortical_separation import GM_LABELS, WM_LABELS


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_legacy_vtk(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read the ASCII POINTS and triangular POLYGONS from a legacy VTK file."""
    lines = path.read_text(encoding="ascii").splitlines()
    point_line = next(i for i, line in enumerate(lines) if line.startswith("POINTS "))
    n_points = int(lines[point_line].split()[1])
    points = np.fromstring(" ".join(lines[point_line + 1:]), sep=" ", count=3 * n_points)
    if points.size != 3 * n_points:
        raise ValueError(f"Could not read all VTK points from {path}")
    polygon_line = next(i for i, line in enumerate(lines) if line.startswith("POLYGONS "))
    n_faces = int(lines[polygon_line].split()[1])
    faces = np.empty((n_faces, 3), dtype=np.int64)
    for index, line in enumerate(lines[polygon_line + 1:polygon_line + 1 + n_faces]):
        values = np.fromstring(line, sep=" ", dtype=np.int64)
        if values.size != 4 or values[0] != 3:
            raise ValueError(f"Expected triangular VTK polygon at {path}:{polygon_line + index + 2}")
        faces[index] = values[1:]
    return points.reshape((-1, 3)), faces


def mesh_topology(path: Path) -> tuple[dict[str, int | bool], np.ndarray, np.ndarray, int]:
    vertices, faces = read_legacy_vtk(path)
    edges = np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]))
    edges.sort(axis=1)
    unique_edges, edge_uses = np.unique(edges, axis=0, return_counts=True)
    adjacency = coo_matrix(
        (np.ones(len(unique_edges), dtype=np.uint8),
         (unique_edges[:, 0], unique_edges[:, 1])),
        shape=(len(vertices), len(vertices)),
    )
    components, component_id = connected_components(adjacency + adjacency.T, directed=False)
    vertex_component = np.bincount(component_id, minlength=components)
    edge_component = component_id[unique_edges[:, 0]]
    face_component = component_id[faces[:, 0]]
    edge_counts_by_component = np.bincount(edge_component, minlength=components)
    face_counts_by_component = np.bincount(face_component, minlength=components)
    largest = int(np.argmax(vertex_component))
    metrics = {
        "vertices": int(len(vertices)),
        "faces": int(len(faces)),
        "edges": int(len(unique_edges)),
        "connected_components": int(components),
        "largest_component_vertices": int(vertex_component[largest]),
        "largest_component_faces": int(face_counts_by_component[largest]),
        "largest_component_euler_characteristic": int(
            vertex_component[largest] - edge_counts_by_component[largest] + face_counts_by_component[largest]
        ),
        "components_under_100_vertices": int(np.count_nonzero(vertex_component < 100)),
        "boundary_edges": int(np.count_nonzero(edge_uses == 1)),
        "nonmanifold_edges": int(np.count_nonzero(edge_uses > 2)),
        "euler_characteristic": int(len(vertices) - len(unique_edges) + len(faces)),
        "watertight_by_edge_incidence": bool(np.all(edge_uses == 2)),
    }
    return metrics, vertices, component_id, largest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("core_dir", type=Path)
    parser.add_argument("--maseg", type=Path, required=True)
    parser.add_argument("--brain-mask", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--watchdog-log", type=Path)
    parser.add_argument("--comparison-core-dir", type=Path)
    args = parser.parse_args()
    core = args.core_dir.resolve()
    surface_dir = core / "surface-thickness"
    files = sorted(
        path for path in core.rglob("*")
        if path.is_file() and path.suffix.lower() in {".gz", ".vtk", ".json"}
    )
    inventory = []
    bad_images = []
    image_cache = {}
    for path in files:
        entry: dict[str, object] = {
            "path": path.relative_to(core).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        if path.name.endswith(".nii.gz"):
            image = nib.load(str(path))
            data = np.asarray(image.dataobj)
            axis_codes = nib.aff2axcodes(image.affine)
            entry.update({
                "shape": list(image.shape),
                "dtype": str(data.dtype),
                "axis_codes": list(axis_codes),
                "finite": bool(np.isfinite(data).all()),
            })
            image_cache[path.stem.removesuffix(".nii")] = image
            if not entry["finite"] or axis_codes != ("R", "A", "S"):
                bad_images.append(entry["path"])
        inventory.append(entry)

    required = {"inner", "central", "outer", "SEG_refine"}
    if not required.issubset(image_cache):
        raise ValueError(f"Missing required core NIfTIs: {sorted(required - set(image_cache))}")
    images = {name: image_cache[name] for name in ("inner", "central", "outer", "SEG_refine")}
    maseg = nib.load(str(args.maseg))
    brain_mask = nib.load(str(args.brain_mask))
    reference = images["inner"]
    for name, image in {**images, "maseg": maseg, "brain_mask": brain_mask}.items():
        if image.shape != reference.shape or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-4):
            raise ValueError(f"{name} does not share the inner level-set grid")
        if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise ValueError(f"{name} is not RAS oriented")
    labels = np.asarray(maseg.dataobj)
    inner = np.asarray(images["inner"].dataobj) <= 0
    central = np.asarray(images["central"].dataobj) <= 0
    outer = np.asarray(images["outer"].dataobj) <= 0
    gm = np.isin(labels, tuple(GM_LABELS))
    wm = np.isin(labels, tuple(WM_LABELS))
    gm_count = int(gm.sum())
    wm_count = int(wm.sum())
    nesting = {
        "inner_outside_central_voxels": int(np.count_nonzero(inner & ~central)),
        "central_outside_outer_voxels": int(np.count_nonzero(central & ~outer)),
    }
    support = {
        "gm_voxels": gm_count,
        "gm_between_inner_outer": int(np.count_nonzero(gm & ~inner & outer)),
        "gm_between_inner_outer_fraction": float(np.count_nonzero(gm & ~inner & outer) / max(1, gm_count)),
        "wm_voxels": wm_count,
        "wm_inside_inner": int(np.count_nonzero(wm & inner)),
        "wm_inside_inner_fraction": float(np.count_nonzero(wm & inner) / max(1, wm_count)),
        "inner_outside_brain_mask_voxels": int(np.count_nonzero(inner & ~(np.asarray(brain_mask.dataobj) > 0))),
    }
    label_support: dict[str, dict[str, int | float | None]] = {}
    for label in sorted(GM_LABELS | WM_LABELS):
        label_mask = labels == label
        voxel_count = int(np.count_nonzero(label_mask))
        if voxel_count == 0:
            continue
        if label in GM_LABELS:
            supported_count = int(np.count_nonzero(label_mask & ~inner & outer))
            region = "gm_between_inner_outer"
        else:
            supported_count = int(np.count_nonzero(label_mask & inner))
            region = "wm_inside_inner"
        label_support[str(label)] = {
            "region": region,
            "voxels": voxel_count,
            "supported_voxels": supported_count,
            "support_fraction": supported_count / voxel_count,
        }
    mesh_results = {
        name: mesh_topology(surface_dir / f"{name}_surface.vtk")
        for name in ("inner", "central", "outer")
    }
    meshes = {name: result[0] for name, result in mesh_results.items()}
    inner_vertices, inner_component_id, largest_inner_component = (
        mesh_results["inner"][1], mesh_results["inner"][2], mesh_results["inner"][3]
    )
    non_dominant_vertices = inner_component_id != largest_inner_component
    vertex_voxels = nib.affines.apply_affine(np.linalg.inv(reference.affine), inner_vertices)
    nearest_voxels = np.rint(vertex_voxels).astype(np.int64)
    inside_grid = np.all((nearest_voxels >= 0) & (nearest_voxels < np.asarray(reference.shape)), axis=1)
    nearest_labels = np.zeros(len(inner_vertices), dtype=labels.dtype)
    nearest_mask = np.zeros(len(inner_vertices), dtype=bool)
    selected_voxels = nearest_voxels[inside_grid]
    nearest_labels[inside_grid] = labels[tuple(selected_voxels.T)]
    nearest_mask[inside_grid] = np.asarray(brain_mask.dataobj)[tuple(selected_voxels.T)] > 0
    non_dominant_in_grid = non_dominant_vertices & inside_grid
    nearest_label_values, nearest_label_counts = np.unique(
        nearest_labels[non_dominant_in_grid], return_counts=True
    )
    shell_locations = {
        "non_dominant_component_count": int(meshes["inner"]["connected_components"] - 1),
        "non_dominant_vertices": int(np.count_nonzero(non_dominant_vertices)),
        "non_dominant_vertices_inside_grid": int(np.count_nonzero(non_dominant_in_grid)),
        "non_dominant_vertices_inside_brain_mask": int(
            np.count_nonzero(non_dominant_in_grid & nearest_mask)
        ),
        "nearest_maseg_labels_inside_grid": {
            str(int(value)): int(count)
            for value, count in zip(nearest_label_values, nearest_label_counts, strict=True)
        },
    }
    thickness = json.loads((surface_dir / "thickness_summary.json").read_text(encoding="utf-8"))
    comparison = None
    if args.comparison_core_dir is not None:
        comparison_root = args.comparison_core_dir.resolve()
        prior_inner_image = nib.load(str(comparison_root / "inner.nii.gz"))
        prior_outer_image = nib.load(str(comparison_root / "outer.nii.gz"))
        for name, image in (("prior inner", prior_inner_image), ("prior outer", prior_outer_image)):
            if image.shape != reference.shape or not np.allclose(
                image.affine, reference.affine, rtol=0, atol=1e-4
            ):
                raise ValueError(f"{name} does not share the current inner level-set grid")
        prior_inner = np.asarray(prior_inner_image.dataobj) <= 0
        prior_outer = np.asarray(prior_outer_image.dataobj) <= 0
        prior_gm = int(np.count_nonzero(gm))
        prior_wm = int(np.count_nonzero(wm))
        prior_gm_support = int(np.count_nonzero(gm & ~prior_inner & prior_outer)) / max(1, prior_gm)
        prior_wm_support = int(np.count_nonzero(wm & prior_inner)) / max(1, prior_wm)
        deltas = []
        for label_text, current in label_support.items():
            label = int(label_text)
            mask = labels == label
            if current["region"] == "gm_between_inner_outer":
                prior_fraction = float(np.count_nonzero(mask & ~prior_inner & prior_outer) / current["voxels"])
            else:
                prior_fraction = float(np.count_nonzero(mask & prior_inner) / current["voxels"])
            deltas.append({
                "label": label,
                "region": current["region"],
                "prior_support_fraction": prior_fraction,
                "current_support_fraction": current["support_fraction"],
                "delta_percentage_points": 100.0 * (current["support_fraction"] - prior_fraction),
            })
        comparison = {
            "core_dir": str(comparison_root),
            "comparison_type": "internal profile comparison, not a source-reference comparison",
            "prior_gm_between_inner_outer_fraction": prior_gm_support,
            "prior_wm_inside_inner_fraction": prior_wm_support,
            "current_minus_prior_gm_percentage_points": 100.0 * (
                support["gm_between_inner_outer_fraction"] - prior_gm_support
            ),
            "current_minus_prior_wm_percentage_points": 100.0 * (
                support["wm_inside_inner_fraction"] - prior_wm_support
            ),
            "largest_per_label_losses": sorted(
                deltas, key=lambda item: item["delta_percentage_points"]
            )[:10],
            "largest_per_label_gains": sorted(
                deltas, key=lambda item: item["delta_percentage_points"], reverse=True
            )[:10],
        }
    log_metrics = None
    if args.watchdog_log is not None:
        log = args.watchdog_log.read_text(encoding="utf-8")
        peaks = [float(value) for value in re.findall(r"peak_sampled_rss_gb=([0-9.]+)", log)]
        log_metrics = {"peak_sampled_rss_gb": max(peaks) if peaks else None}
    report = {
        "evidence_category": "in-vivo integration and quality audit; not source equivalence or clinical validation",
        "core_dir": str(core),
        "input_maseg": str(args.maseg.resolve()),
        "input_brain_mask": str(args.brain_mask.resolve()),
        "n_files_hashed": len(inventory),
        "n_nifti": sum(path.name.endswith(".nii.gz") for path in files),
        "bad_nifti": bad_images,
        "files": inventory,
        "maseg_labels_208_209": {
            str(value): int(np.count_nonzero(labels == value)) for value in (208, 209)
        },
        "quality": {**support, **nesting},
        "label_support_by_label": label_support,
        "comparison": comparison,
        "meshes": meshes,
        "inner_non_dominant_shell_locations": shell_locations,
        "thickness": thickness,
        "watchdog": log_metrics,
        "source_reference_comparison": "unavailable: supplied output directory had no historical reference products",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {args.report}; {len(inventory)} files hashed; bad_nifti={len(bad_images)}")


if __name__ == "__main__":
    main()
