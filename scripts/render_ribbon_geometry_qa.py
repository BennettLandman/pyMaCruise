#!/usr/bin/env python3
"""Render a QA-only multiplanar overlay and 3D surface view for a base run.

This script is for visual inspection, not a production visualization feature or
quantitative acceptance test. It reads completed run products and writes PNGs
and a small provenance JSON to an ignored output directory.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import nibabel as nib
import numpy as np
from scipy.ndimage import map_coordinates


GM_LABELS = np.array([
    *range(100, 110), *range(112, 126), 128, 129,
    *range(132, 158), *range(160, 188), *range(190, 208),
], dtype=np.uint16)
WM_LABELS = np.array([4, 23, 30, 36, 37, 44, 45, 51, 52, *range(55, 63)], dtype=np.uint16)
SURFACE_NAMES = ("inner", "central", "outer")
SURFACE_COLORS = {"inner": "#ff453a", "central": "#ffd60a", "outer": "#30d158"}


def _outputs_dir(manifest: dict, stage: str) -> Path:
    outputs = manifest["stages"][stage]["outputs"]
    paths = [Path(item["path"]) for item in outputs if "path" in item]
    if not paths:
        raise ValueError(f"Manifest stage {stage!r} has no file outputs")
    # All files for a stage live in its attempt directory.
    return paths[0].parent


def resolve_products(manifest_or_attempt: Path) -> tuple[Path, Path]:
    if manifest_or_attempt.is_file():
        manifest = json.loads(manifest_or_attempt.read_text(encoding="utf-8"))
        return _outputs_dir(manifest, "cruise-core"), _outputs_dir(manifest, "surface-thickness")
    if manifest_or_attempt.is_dir():
        root = manifest_or_attempt
        if (root / "run.json").is_file():
            return resolve_products(root / "run.json")
        # A stage attempt directory can be either the core or surface directory.
        core_candidates = list(root.rglob("inner.nii.gz"))
        surface_candidates = list(root.rglob("inner_surface.vtk"))
        if core_candidates and surface_candidates:
            return core_candidates[0].parent, surface_candidates[0].parent
        parent = root.parent.parent if root.name.startswith("attempt-") else root
        core_candidates = list(parent.rglob("inner.nii.gz"))
        surface_candidates = list(parent.rglob("inner_surface.vtk"))
        if core_candidates and surface_candidates:
            return core_candidates[-1].parent, surface_candidates[-1].parent
    raise ValueError("Supply a run.json manifest or a directory containing the completed base-run attempts")


def read_mesh(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a legacy ASCII VTK triangle mesh."""
    points: np.ndarray | None = None
    faces: list[tuple[int, int, int]] = []
    with path.open("r", encoding="ascii") as stream:
        lines = iter(stream)
        for line in lines:
            if line.startswith("POINTS "):
                n_points = int(line.split()[1])
                values: list[float] = []
                while len(values) < 3 * n_points:
                    values.extend(float(v) for v in next(lines).split())
                points = np.asarray(values, dtype=np.float32).reshape((-1, 3))
            elif line.startswith("POLYGONS "):
                n_faces = int(line.split()[1])
                for idx in range(n_faces):
                    row = next(lines).split()
                    if len(row) != 4 or row[0] != "3":
                        raise ValueError(f"Expected triangular VTK polygons in {path}")
                    faces.append((int(row[1]), int(row[2]), int(row[3])))
                if points is None:
                    raise ValueError(f"VTK has polygons before points: {path}")
                return points, np.asarray(faces, dtype=np.int32)
    raise ValueError(f"No triangular POLYGONS found in {path}")


def _sample_plane(image: nib.spatialimages.SpatialImage, reference: nib.spatialimages.SpatialImage,
                  axis: int, index: int, order: int, dtype: np.dtype) -> np.ndarray:
    """Sample an input volume on an orthogonal plane of the output level-set grid."""
    other = [i for i in range(3) if i != axis]
    shape = reference.shape[:3]
    u = np.arange(shape[other[0]], dtype=np.float32)
    v = np.arange(shape[other[1]], dtype=np.float32)
    uu, vv = np.meshgrid(u, v, indexing="ij")
    voxel = np.zeros((3, uu.size), dtype=np.float64)
    voxel[axis] = index
    voxel[other[0]] = uu.ravel()
    voxel[other[1]] = vv.ravel()
    world = reference.affine @ np.vstack((voxel, np.ones(uu.size)))
    source = np.linalg.inv(image.affine) @ world
    # Load one proxy-backed array once for each plane. Human inputs are modest
    # (typically tens of megabytes); only the selected 2D plane is interpolated.
    volume = np.asanyarray(image.dataobj)
    sampled = map_coordinates(volume, source[:3], order=order, mode="constant", cval=0, prefilter=False)
    return sampled.reshape(uu.shape).astype(dtype, copy=False)


def render_slices(core: Path, t1_path: Path, maseg_path: Path, out: Path,
                  indices: tuple[int, int, int], title: str) -> dict:
    reference = nib.load(str(core / "inner.nii.gz"))
    t1 = nib.load(str(t1_path))
    maseg = nib.load(str(maseg_path))
    if len(reference.shape) != 3 or len(t1.shape) != 3 or len(maseg.shape) != 3:
        raise ValueError("T1, MaSeg, and level-set images must be 3D")
    ras = nib.aff2axcodes(reference.affine)
    if ras != ("R", "A", "S"):
        raise ValueError(f"Surface level-set grid is {ras}, expected RAS")
    for name in SURFACE_NAMES:
        product = nib.load(str(core / f"{name}.nii.gz"))
        if product.shape != reference.shape or not np.allclose(product.affine, reference.affine):
            raise ValueError(f"{name} level set does not share inner grid")

    fig, axes = plt.subplots(1, 3, figsize=(17, 6), constrained_layout=True)
    descriptions = ((0, indices[0], "Sagittal (R index)"),
                    (1, indices[1], "Coronal (A index)"),
                    (2, indices[2], "Axial (S index)"))
    levelsets = [nib.load(str(core / f"{name}.nii.gz")) for name in SURFACE_NAMES]
    for ax, (axis, index, label) in zip(axes, descriptions):
        image_plane = _sample_plane(t1, reference, axis, index, 1, np.float32)
        maseg_plane = _sample_plane(maseg, reference, axis, index, 0, np.uint16)
        low, high = np.percentile(image_plane[image_plane > 0], [1, 99]) if np.any(image_plane > 0) else (0, 1)
        ax.imshow(image_plane.T, cmap="gray", origin="lower", vmin=low, vmax=high)
        gm = np.isin(maseg_plane, GM_LABELS)
        wm = np.isin(maseg_plane, WM_LABELS)
        if gm.any():
            ax.contour(gm.T.astype(np.uint8), levels=[0.5], colors=["#40c8e0"], linewidths=0.55, alpha=0.8)
        if wm.any():
            ax.contour(wm.T.astype(np.uint8), levels=[0.5], colors=["#bf5af2"], linewidths=0.55, alpha=0.8)
        for name, image in zip(SURFACE_NAMES, levelsets):
            plane = _sample_plane(image, reference, axis, index, 1, np.float32)
            ax.contour(plane.T, levels=[0], colors=[SURFACE_COLORS[name]], linewidths=1.0)
        ax.set_title(f"{label}: {index}")
        ax.set_axis_off()
    fig.suptitle(title + "\nRAS; surfaces: inner red, central yellow, outer green; MaSeg GM cyan / WM purple", fontsize=11)
    path = out / "ribbon_contours_multiplanar.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return {"file": path.name, "slice_indices_ras_voxel": {"sagittal_r": indices[0], "coronal_a": indices[1], "axial_s": indices[2]}}


def render_meshes(surface_dir: Path, out: Path, title: str, spacing_mm: float) -> dict:
    fig = plt.figure(figsize=(16, 6), constrained_layout=True)
    mesh_meta = {}
    bounds = []
    loaded = {}
    for name in SURFACE_NAMES:
        vertices, faces = read_mesh(surface_dir / f"{name}_surface.vtk")
        bounds.append((vertices.min(axis=0), vertices.max(axis=0)))
        # Cluster nearby vertices onto a display grid, then collapse duplicate
        # and degenerate triangles. Unlike uniform face striding, this retains
        # a connected-looking coarse skin on dense meshes.
        quantized = np.rint(vertices / spacing_mm).astype(np.int32)
        quantized_unique, inverse = np.unique(quantized, axis=0, return_inverse=True)
        display_vertices = quantized_unique.astype(np.float32) * spacing_mm
        display_faces = inverse[faces]
        valid = ((display_faces[:, 0] != display_faces[:, 1]) &
                 (display_faces[:, 1] != display_faces[:, 2]) &
                 (display_faces[:, 2] != display_faces[:, 0]))
        display_faces = display_faces[valid]
        display_faces.sort(axis=1)
        display_faces = np.unique(display_faces, axis=0)
        loaded[name] = (display_vertices, display_faces)
        mesh_meta[name] = {"source_faces": int(len(faces)), "display_faces": int(len(display_faces)),
                           "source_vertices": int(len(vertices)), "display_vertices": int(len(display_vertices)),
                           "display_vertex_cluster_spacing_mm": spacing_mm}
        del quantized, quantized_unique, inverse, faces, vertices
    low = np.min([b[0] for b in bounds], axis=0)
    high = np.max([b[1] for b in bounds], axis=0)
    center = (low + high) / 2
    radius = float(np.max(high - low) / 2)
    for i, name in enumerate(SURFACE_NAMES, 1):
        ax = fig.add_subplot(1, 3, i, projection="3d")
        vertices, faces = loaded[name]
        ax.add_collection3d(Poly3DCollection(vertices[faces], facecolor=SURFACE_COLORS[name], edgecolor="none", alpha=0.92))
        ax.set_xlim(center[0] - radius, center[0] + radius)
        ax.set_ylim(center[1] - radius, center[1] + radius)
        ax.set_zlim(center[2] - radius, center[2] + radius)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("R (mm)")
        ax.set_ylabel("A (mm)")
        ax.set_zlabel("S (mm)")
        ax.set_title(name)
        del vertices, faces
    fig.suptitle(title + "\n3D meshes shown separately; display-only vertex clustering", fontsize=12)
    path = out / "ribbon_surfaces_3d.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return {"file": path.name, "meshes": mesh_meta, "display_only_vertex_clustering": True}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_or_attempt", type=Path, help="base run.json manifest, run folder, or stage attempt directory")
    parser.add_argument("--t1", type=Path, required=True, help="RAS T1 NIfTI (.nii.gz preferred)")
    parser.add_argument("--maseg", type=Path, required=True, help="RAS MaSeg NIfTI (.nii.gz preferred)")
    parser.add_argument("--output-dir", type=Path, required=True, help="Ignored QA output folder, preferably under build/")
    parser.add_argument("--slices", type=int, nargs=3, metavar=("R", "A", "S"), help="Optional sagittal, coronal, axial indices on the run grid")
    parser.add_argument("--label", default="HISTORICAL / STALE OUTPUTS", help="Visible title label; default prevents stale outputs being mistaken for current")
    parser.add_argument("--display-spacing-mm", type=float, default=1.25, help="Vertex clustering spacing for 3D display only")
    args = parser.parse_args()
    core, surface = resolve_products(args.run_or_attempt.resolve())
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    reference = nib.load(str(core / "inner.nii.gz"))
    slices = args.slices or tuple(int((n - 1) // 2) for n in reference.shape[:3])
    title = args.label
    slice_result = render_slices(core, args.t1.resolve(), args.maseg.resolve(), out, slices, title)
    mesh_result = render_meshes(surface, out, title, args.display_spacing_mm)
    provenance = {
        "label": args.label,
        "status_warning": "Visualization QA only; inspect the corresponding run manifest before interpreting as current.",
        "core_attempt_dir": str(core), "surface_attempt_dir": str(surface),
        "t1": str(args.t1.resolve()), "maseg": str(args.maseg.resolve()),
        "ras_axis_codes": list(nib.aff2axcodes(reference.affine)),
        "slice_view": slice_result, "mesh_view": mesh_result,
    }
    (out / "qa_render_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
