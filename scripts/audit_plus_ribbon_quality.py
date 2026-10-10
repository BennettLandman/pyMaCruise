#!/usr/bin/env python3
"""Measure ribbon quality on outputs from a completed plus workflow manifest.

This audit deliberately reports failed quality criteria instead of treating a
completed workflow manifest as evidence that its synthetic anatomy is useful.
"""

from __future__ import annotations

import argparse
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

from pymacruise.surfaces import extract_levelset_surface  # noqa: E402
from pymacruise.thickness import closest_point_thickness, levelset_subtraction_thickness  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="completed plus-full-workflow.json")
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        manifest.get("status") != "complete"
        or manifest.get("terminal") is not True
        or manifest.get("stage") not in {
            "plus-full-workflow",
            "synthetic-plus-terminal-surface-composition",
        }
    ):
        raise ValueError("ribbon audit requires a completed terminal plus or terminal-surface-composition manifest")
    outputs = manifest["final_outputs"]
    image_names = ("plus_inner_level_set", "plus_central_level_set", "plus_outer_level_set")
    images = {}
    for key in image_names:
        record = outputs.get(key, {})
        image_path = Path(record.get("path", "")).resolve()
        if (
            not image_path.is_file()
            or image_path.is_symlink()
            or hashlib.sha256(image_path.read_bytes()).hexdigest() != record.get("sha256")
        ):
            raise ValueError(f"final output {key!r} is missing or fails its manifest hash")
        images[key] = nib.load(str(image_path))
    labels_path = manifest["inputs"]["maseg"]["path"]
    if (
        not Path(labels_path).is_file()
        or hashlib.sha256(Path(labels_path).read_bytes()).hexdigest()
        != manifest["inputs"]["maseg"].get("sha256")
    ):
        raise ValueError("MaSeg input is missing or fails its manifest hash")
    labels_image = nib.load(labels_path)
    ref = images["plus_inner_level_set"]
    for key, image in {**images, "maseg": labels_image}.items():
        if image.shape != ref.shape or not np.allclose(image.affine, ref.affine):
            raise ValueError(f"{key} does not share the final surface grid")
        if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise ValueError(f"{key} is not on a RAS grid")
        if not np.isfinite(np.asarray(image.dataobj)).all():
            raise ValueError(f"{key} contains non-finite values")

    surfaces = {
        key: extract_levelset_surface(image)
        for key, image in zip(("inner", "central", "outer"), images.values(), strict=True)
    }
    topology = {}
    for name, mesh in surfaces.items():
        edges = np.concatenate(
            (mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]])
        )
        edges.sort(axis=1)
        unique_edges, edge_uses = np.unique(edges, axis=0, return_counts=True)
        adjacency = coo_matrix(
            (np.ones(len(unique_edges), dtype=np.uint8),
             (unique_edges[:, 0], unique_edges[:, 1])),
            shape=(len(mesh.vertices), len(mesh.vertices)),
        )
        component_count, _ = connected_components(adjacency + adjacency.T, directed=False)
        topology[name] = {
            "connected_components": int(component_count),
            "boundary_edges": int(np.count_nonzero(edge_uses == 1)),
            "nonmanifold_edges": int(np.count_nonzero(edge_uses > 2)),
            "euler_characteristic": int(len(mesh.vertices) - len(unique_edges) + len(mesh.faces)),
            "watertight": bool(np.all(edge_uses == 2)),
        }
    thickness = closest_point_thickness(
        surfaces["inner"], surfaces["outer"], central=surfaces["central"]
    )
    vertex_values = thickness.inner_surface.point_data["mean_thickness_mm"]
    volume = levelset_subtraction_thickness(
        images["plus_inner_level_set"], images["plus_outer_level_set"]
    )
    label_values, label_counts = np.unique(
        np.asarray(labels_image.dataobj), return_counts=True
    )
    expected_labels = {0, 5, 44, 100}
    labels = np.asarray(labels_image.dataobj)
    inner_region = np.asarray(images["plus_inner_level_set"].dataobj) <= 0
    central_region = np.asarray(images["plus_central_level_set"].dataobj) <= 0
    outer_region = np.asarray(images["plus_outer_level_set"].dataobj) <= 0
    nesting = {
        "inner_outside_central_voxels": int(np.count_nonzero(inner_region & ~central_region)),
        "central_outside_outer_voxels": int(np.count_nonzero(central_region & ~outer_region)),
    }
    wm = labels == 44
    gm = labels == 100
    csf = labels == 5
    label_disagreements = {
        "wm_voxels_outside_inner": int(np.count_nonzero(wm & ~inner_region)),
        "gm_voxels_inside_inner": int(np.count_nonzero(gm & inner_region)),
        "gm_voxels_outside_outer": int(np.count_nonzero(gm & ~outer_region)),
        "csf_voxels_inside_outer": int(np.count_nonzero(csf & outer_region)),
    }
    label_support = {
        "wm_inside_inner": 1.0 - label_disagreements["wm_voxels_outside_inner"] / max(1, int(np.count_nonzero(wm))),
        "gm_between_inner_and_outer": (
            int(np.count_nonzero(gm))
            - label_disagreements["gm_voxels_inside_inner"]
            - label_disagreements["gm_voxels_outside_outer"]
        ) / max(1, int(np.count_nonzero(gm))),
        "csf_outside_outer": 1.0 - label_disagreements["csf_voxels_inside_outer"] / max(1, int(np.count_nonzero(csf))),
    }
    quality = {
        "nominal_ribbon_mm": 5.0,
        "minimum_label_support": 0.97,
        "surface_median_tolerance_mm": 1.0,
        "surface_p05_minimum_mm": 2.0,
        "surface_p95_maximum_mm": 8.0,
        "levelset_mean_tolerance_mm": 1.0,
    }
    surface_median = float(np.median(vertex_values))
    surface_p05, surface_p95 = map(float, np.percentile(vertex_values, [5, 95]))
    label_set = {int(value) for value in label_values}
    checks = {
        "surfaces_closed_connected_euler_two": all(
            metric["watertight"] and metric["boundary_edges"] == 0
            and metric["nonmanifold_edges"] == 0
            and metric["connected_components"] == 1
            and metric["euler_characteristic"] == 2
            for metric in topology.values()
        ),
        "explicit_maseg_wm_gm_csf_labels_present": expected_labels.issubset(label_set),
        "labeled_tissue_regions_agree_with_surfaces": bool(
            min(label_support.values()) >= quality["minimum_label_support"]
        ),
        "nested_voxel_regions": all(value == 0 for value in nesting.values()),
        "production_surface_thickness_within_bounds": bool(
            abs(surface_median - quality["nominal_ribbon_mm"])
            <= quality["surface_median_tolerance_mm"]
            and surface_p05 >= quality["surface_p05_minimum_mm"]
            and surface_p95 <= quality["surface_p95_maximum_mm"]
        ),
        "levelset_mean_thickness_within_bounds": bool(
            volume.crossing_voxel_count > 0
            and np.isfinite(volume.mean_mm)
            and abs(volume.mean_mm - quality["nominal_ribbon_mm"])
            <= quality["levelset_mean_tolerance_mm"]
        ),
    }
    report = {
        "evidence_category": "synthetic plus workflow quality audit; not Java equivalence or real-data validation",
        "manifest": str(manifest_path.relative_to(ROOT))
        if manifest_path.is_relative_to(ROOT) else str(manifest_path),
        "manifest_run_signature": manifest["run_signature"],
        "shape": list(ref.shape),
        "orientation": "RAS",
        "maseg_labels": {
            "present": [int(value) for value in label_values],
            "counts": {str(int(value)): int(count) for value, count in zip(label_values, label_counts, strict=True)},
        },
        "label_surface_consistency": {
            "disagreement_voxel_counts": label_disagreements,
            "support_fractions": label_support,
            "minimum_required_support": quality["minimum_label_support"],
        },
        "nesting_voxel_counts": nesting,
        "surface_meshes": {
            name: {
                "vertices": int(len(mesh.vertices)),
                "faces": int(len(mesh.faces)),
                **topology[name],
            }
            for name, mesh in surfaces.items()
        },
        "surface_closest_point_thickness": {
            "mean_mm": float(thickness.mean_mm),
            "stdev_mm": float(thickness.standard_deviation_mm),
            "median_mm": surface_median,
            "p05_mm": surface_p05,
            "p95_mm": surface_p95,
            "sample_count": int(thickness.nonzero_vertex_count),
        },
        "levelset_subtraction_thickness": {
            "crossing_voxel_count": int(volume.crossing_voxel_count),
            "mean_mm": float(volume.mean_mm),
            "stdev_mm": float(volume.standard_deviation_mm),
            "minimum_mm": float(volume.minimum_mm),
            "maximum_mm": float(volume.maximum_mm),
        },
        "quality_bounds": quality,
        "checks": checks,
        "passed": all(checks.values()),
    }
    report_path = args.report.resolve()
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
