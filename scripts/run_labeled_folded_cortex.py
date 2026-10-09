#!/usr/bin/env python3
"""Run source-selected MaCRUISE TGDM on the labeled synthetic cortex.

This is synthetic integration evidence. It does not claim Java equivalence or
real-image validation. All inputs and generated level sets are RAS NIfTIs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pymacruise.tgdm import TGDMBaselineParameters, tgdm_nested_arrays  # noqa: E402
from pymacruise.surfaces import SurfaceMesh, extract_levelset_surface  # noqa: E402
from pymacruise.thickness import (  # noqa: E402
    closest_point_thickness,
    levelset_subtraction_thickness,
)


def _mesh_metrics(surface: SurfaceMesh) -> dict[str, object]:
    vertices, faces = surface.vertices, surface.faces
    edges = np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]))
    edges.sort(axis=1)
    unique_edges, edge_uses = np.unique(edges, axis=0, return_counts=True)
    adjacency = coo_matrix(
        (np.ones(len(unique_edges), dtype=np.uint8),
         (unique_edges[:, 0], unique_edges[:, 1])),
        shape=(len(vertices), len(vertices)),
    )
    adjacency = adjacency + adjacency.T
    components, _ = connected_components(adjacency, directed=False)
    euler = len(vertices) - len(unique_edges) + len(faces)
    return {
        "vertices": int(len(vertices)),
        "faces": int(len(faces)),
        "connected_components": int(components),
        "boundary_edges": int(np.count_nonzero(edge_uses == 1)),
        "nonmanifold_edges": int(np.count_nonzero(edge_uses > 2)),
        "euler_characteristic": int(euler),
        "watertight": bool(np.all(edge_uses == 2)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fixture", type=Path,
        default=ROOT / "tests/data/synthetic/labeled_folded_cortex",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    names = (
        "initial_membership", "white_matter", "gray_matter", "ace_gray_matter",
        "ace_skeleton", "brain_mask", "maseg_labels", "t1",
    )
    images = {name: nib.load(str(args.fixture / f"{name}.nii.gz")) for name in names}
    reference = images["initial_membership"]
    if nib.aff2axcodes(reference.affine) != ("R", "A", "S"):
        raise ValueError("fixture is not RAS")
    for name, image in images.items():
        if image.shape != reference.shape or not np.allclose(image.affine, reference.affine):
            raise ValueError(f"{name} does not share the initial-membership grid")
        if not np.isfinite(np.asarray(image.dataobj)).all():
            raise ValueError(f"{name} contains non-finite values")

    data = {name: np.asarray(image.dataobj) for name, image in images.items()}
    labels = data["maseg_labels"]
    label_values = set(int(value) for value in np.unique(labels))
    if not {0, 5, 44, 100}.issubset(label_values):
        raise ValueError(f"expected explicit MaSeg labels 0/5/44/100, got {sorted(label_values)}")

    # This is the complete selected source profile. `run_ToadsCruise` changes
    # the wrapper default WM iso from 0.5 to 0.6; all iteration counts remain
    # at their selected source defaults. The selected Java wrapper supplies a
    # zero GVF to TGDM's central stage.
    params = TGDMBaselineParameters(white_matter_iso=0.6)
    result = tgdm_nested_arrays(
        data["initial_membership"], data["white_matter"], data["gray_matter"],
        data["ace_gray_matter"], brain_mask=data["brain_mask"],
        ace_skeleton=data["ace_skeleton"],
        gvf=np.zeros((*reference.shape, 3), dtype=np.float32),
        parameters=params,
    )

    fields = {
        "inner": result.inner.level_set,
        "central": result.central.level_set,
        "outer": result.outer.level_set,
    }
    all_fields = {
        "signed_seed": result.initial_level_set,
        "smoothed_initial": result.smoothed_initial.level_set,
        "initial_inner": result.initial_inner.level_set,
        "tissue_sign": result.tissue_sign,
        "central": result.central.level_set,
        "central_shrink": result.central_shrink.level_set,
        "inner": result.inner.level_set,
        "outer": result.outer.level_set,
    }
    for name, field in all_fields.items():
        out = nib.Nifti1Image(np.asarray(field), reference.affine, reference.header.copy())
        out.set_qform(reference.affine, code=1)
        out.set_sform(reference.affine, code=1)
        # The initial-membership fixture is uint8; keep the NIfTI header's
        # stored datatype aligned to the stage array or float level sets are
        # silently quantized when written through its copied header.
        out.set_data_dtype(np.asarray(field).dtype)
        output_path = args.output_dir / f"{name}.nii.gz"
        nib.save(out, str(output_path))
        if nib.load(str(output_path)).get_data_dtype() != np.asarray(field).dtype:
            raise ValueError(f"{name} output NIfTI datatype does not match its array")
    meshes = {}
    surfaces = {}
    levelset_images = {}
    for name, field in fields.items():
        if not np.isfinite(field).all():
            raise ValueError(f"{name} level set contains non-finite values")
        out = nib.Nifti1Image(field.astype(np.float32), reference.affine, reference.header.copy())
        out.set_qform(reference.affine, code=1)
        out.set_sform(reference.affine, code=1)
        out.set_data_dtype(np.float32)
        output_path = args.output_dir / f"{name}_levelset.nii.gz"
        nib.save(out, str(output_path))
        if nib.load(str(output_path)).get_data_dtype() != np.dtype(np.float32):
            raise ValueError(f"{name} level-set NIfTI was not stored as float32")
        levelset_image = nib.load(str(output_path))
        surface = extract_levelset_surface(levelset_image)
        surfaces[name] = surface
        levelset_images[name] = levelset_image
        meshes[name] = _mesh_metrics(surface)
        meshes[name]["extractor"] = "PyMaCRUISE Java-derived JIST IsoSurfaceOnGrid connectivity-18 templates"
        mesh = meshes[name]
        if mesh["boundary_edges"] != 0 or mesh["nonmanifold_edges"] != 0:
            raise ValueError(f"{name} surface is open or nonmanifold")

    inner_mask = fields["inner"] <= 0
    central_mask = fields["central"] <= 0
    outer_mask = fields["outer"] <= 0
    nesting = {
        "inner_outside_central_voxels": int(np.count_nonzero(inner_mask & ~central_mask)),
        "central_outside_outer_voxels": int(np.count_nonzero(central_mask & ~outer_mask)),
    }
    # Label support checks use signed-distance level-set signs and report the
    # exact disagreements; one-voxel boundary bands are separately counted.
    wm = labels == 44
    gm = labels == 100
    csf = labels == 5
    label_check = {
        "wm_voxels_outside_inner": int(np.count_nonzero(wm & ~inner_mask)),
        "gm_voxels_inside_inner": int(np.count_nonzero(gm & inner_mask)),
        "gm_voxels_outside_outer": int(np.count_nonzero(gm & ~outer_mask)),
        "csf_voxels_inside_outer": int(np.count_nonzero(csf & outer_mask)),
    }
    wm_inside_fraction = 1.0 - label_check["wm_voxels_outside_inner"] / max(1, int(np.count_nonzero(wm)))
    gm_between_count = int(np.count_nonzero(gm)) - label_check["gm_voxels_inside_inner"] - label_check["gm_voxels_outside_outer"]
    gm_between_fraction = gm_between_count / max(1, int(np.count_nonzero(gm)))
    csf_outside_fraction = 1.0 - label_check["csf_voxels_inside_outer"] / max(1, int(np.count_nonzero(csf)))
    # Vertex-to-vertex nearest distances are a reproducible sampled ribbon
    # thickness diagnostic. They are not an exact point-to-triangle distance.
    inner_vertices = surfaces["inner"].vertices
    outer_vertices = surfaces["outer"].vertices
    inner_to_outer = cKDTree(outer_vertices).query(inner_vertices, k=1)[0]
    outer_to_inner = cKDTree(inner_vertices).query(outer_vertices, k=1)[0]
    ribbon_distances = np.concatenate((inner_to_outer, outer_to_inner))
    if not np.isfinite(ribbon_distances).all() or np.any(ribbon_distances <= 0):
        raise ValueError("ribbon surface distance is not finite and strictly positive")

    # Validate the production triangle-distance thickness path as well as the
    # source level-set subtraction output. These are distinct measures and are
    # both retained in the report rather than blended into one tolerance.
    surface_thickness = closest_point_thickness(
        surfaces["inner"], surfaces["outer"], central=surfaces["central"]
    )
    point_thickness = surface_thickness.inner_surface.point_data["mean_thickness_mm"]
    volume_thickness = levelset_subtraction_thickness(
        levelset_images["inner"], levelset_images["outer"]
    )
    thickness_path = args.output_dir / "thickness_levelset.nii.gz"
    nib.save(volume_thickness.image, str(thickness_path))
    if nib.load(str(thickness_path)).get_data_dtype() != np.dtype(np.float32):
        raise ValueError("level-set subtraction thickness output was not stored as float32")
    if (not np.isfinite(point_thickness).all() or np.any(point_thickness <= 0)
            or volume_thickness.crossing_voxel_count == 0
            or not np.isfinite(volume_thickness.mean_mm)):
        raise ValueError("production ribbon thickness contains invalid or empty measurements")

    checks = {
        "all_inputs_share_finite_ras_grid": True,
        "explicit_maseg_labels_present": True,
        "all_surfaces_closed_manifold": all(
            mesh["watertight"] and mesh["boundary_edges"] == 0
            and mesh["nonmanifold_edges"] == 0 for mesh in meshes.values()
        ),
        "each_surface_one_component_euler_2": all(
            mesh["connected_components"] == 1 and mesh["euler_characteristic"] == 2
            for mesh in meshes.values()
        ),
        "nested_voxel_regions": all(value == 0 for value in nesting.values()),
        "positive_finite_inner_outer_vertex_distances": bool(
            np.isfinite(ribbon_distances).all() and np.all(ribbon_distances > 0)
        ),
        "production_surface_and_volume_thickness_valid": bool(
            np.isfinite(point_thickness).all() and np.all(point_thickness > 0)
            and volume_thickness.crossing_voxel_count > 0
            and np.isfinite(volume_thickness.mean_mm)
        ),
        "labeled_tissue_regions_agree_with_surfaces": bool(
            wm_inside_fraction >= 0.97 and gm_between_fraction >= 0.97
            and csf_outside_fraction >= 0.97
        ),
        "nominal_five_mm_ribbon_thickness_retained": bool(
            abs(float(np.median(ribbon_distances)) - 5.0) <= 1.0
            and float(np.percentile(ribbon_distances, 5)) >= 2.0
            and float(np.percentile(ribbon_distances, 95)) <= 8.0
        ),
    }
    report = {
        "evidence_category": "synthetic integration; not Java equivalence or real-data validation",
        "fixture": str(args.fixture.relative_to(ROOT)),
        "shape": list(reference.shape),
        "orientation": "RAS",
        "affine_ras_mm": reference.affine.tolist(),
        "labels": {"present": sorted(label_values), "counts": {
            str(value): int(np.count_nonzero(labels == value)) for value in sorted(label_values)
        }},
        "source_selected_parameters": {
            "initial_iso_level": params.initial_iso_level,
            "white_matter_iso": params.white_matter_iso,
            "gray_matter_iso": params.gray_matter_iso,
            "initial_smoothing": {"curvature": params.initial_curvature, "iterations": params.initial_iterations},
            "inner": {"curvature": params.inner_curvature, "pressure": params.inner_pressure,
                      "external": params.inner_external, "iterations": params.inner_iterations},
            "central": {"curvature": params.central_curvature, "pressure": params.central_pressure,
                        "external": params.central_external, "iterations": params.central_iterations},
            "central_shrink_iterations": params.central_shrink_iterations,
            "outer": {"curvature": params.outer_curvature, "pressure": params.outer_pressure,
                      "external": params.outer_external,
                      "effective_iterations": params.central_iterations,
                      "wrapper_note": "selected wrapper passes central iteration option to Java outer stage"},
            "central_gvf": "zero vector, as selected wrapper receives from absent GVF stage in this controlled fixture",
        },
        "completed_steps": {
            "smoothing": result.smoothed_initial.completed_steps,
            "initial_inner": result.initial_inner.completed_steps,
            "central": result.central.completed_steps,
            "central_shrink": result.central_shrink.completed_steps,
            "inner": result.inner.completed_steps,
            "outer": result.outer.completed_steps,
        },
        "surfaces": meshes,
        "nesting_voxel_counts": nesting,
        "label_surface_consistency_counts": label_check,
        "label_surface_consistency_fractions": {
            "wm_inside_inner": wm_inside_fraction,
            "gm_between_inner_and_outer": gm_between_fraction,
            "csf_outside_outer": csf_outside_fraction,
            "minimum": 0.97,
        },
        "ribbon_thickness_vertex_nearest_mm": {
            "minimum": float(np.min(ribbon_distances)),
            "p05": float(np.percentile(ribbon_distances, 5)),
            "median": float(np.median(ribbon_distances)),
            "p95": float(np.percentile(ribbon_distances, 95)),
            "maximum": float(np.max(ribbon_distances)),
            "sampling_note": "symmetric nearest vertex distance; approximate diagnostic, not point-to-triangle thickness",
        },
        "production_ribbon_thickness": {
            "surface_closest_point": {
                "mean_mm": float(surface_thickness.mean_mm),
                "stdev_mm": float(surface_thickness.standard_deviation_mm),
                "sample_count": int(surface_thickness.nonzero_vertex_count),
                "minimum_mm": float(np.min(point_thickness)),
                "p05_mm": float(np.percentile(point_thickness, 5)),
                "median_mm": float(np.median(point_thickness)),
                "p95_mm": float(np.percentile(point_thickness, 95)),
                "maximum_mm": float(np.max(point_thickness)),
            },
            "levelset_subtraction": {
                "crossing_voxels": int(volume_thickness.crossing_voxel_count),
                "mean_mm": float(volume_thickness.mean_mm),
                "stdev_mm": float(volume_thickness.standard_deviation_mm),
                "minimum_mm": float(volume_thickness.minimum_mm),
                "maximum_mm": float(volume_thickness.maximum_mm),
                "output": str(thickness_path.relative_to(ROOT)),
            },
            "implementation": "production PyMaCRUISE closest-point triangle and level-set subtraction APIs",
        },
        "checks": checks,
        "passed": all(checks.values()),
        "java_acceptance_gates_not_yet_evaluated": {
            "minimum_region_dice": 0.97,
            "maximum_symmetric_mean_surface_distance_mm": 0.01,
            "maximum_symmetric_hausdorff_mm": 1.0,
        },
    }
    (args.output_dir / "ribbon-quality.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if not report["passed"]:
        raise SystemExit("one or more source-selected synthetic ribbon checks failed")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
