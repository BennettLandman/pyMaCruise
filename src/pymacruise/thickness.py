"""Source-mapped CRUISE level-set and closest-point thickness measurements.

The selected wrappers call VandyCRUISE ``LevelSetThickness`` with its only
enabled method (level-set subtraction), and ``ClosestPointThicknessForInner``
for the legacy surface statistic. Surface distances use exact closest points
on triangles, accelerated by a centroid k-d tree with a conservative bound.
Source versions: ``LevelSetThickness.java`` CVS r1.2,
``GetThicknessCRUISE.java`` CVS r1.2, and
``ClosestPointThicknessForInner.java`` CVS r1.2 (author Yuankai Huo).

Paper context: Han et al. (2004), Eq. 13, defines CRUISE thickness at a point
on the triangulated central surface as the difference between the inner and
outer signed-distance level sets there. The selected wrappers port a volume
level-set subtraction product and a separate closest-point surface measure;
neither adapter is claimed to reproduce that central-surface equation or the
paper's measurements. The different JIST-to-Lewiner surface substitution also
prevents an assumption of equivalent mesh-dependent values. See
``docs/papers/cruise.md`` for the operation-to-code and validation boundaries.
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras
from pymacruise.jist_intersector import JistSurfaceTree
from pymacruise.surfaces import SurfaceMesh


@dataclass(frozen=True)
class ThicknessVolumeResult:
    """Level-set subtraction volume and statistics in physical millimetres."""

    image: nib.Nifti1Image
    mean_mm: float
    standard_deviation_mm: float
    minimum_mm: float
    maximum_mm: float
    crossing_voxel_count: int


@dataclass(frozen=True)
class SurfaceThicknessResult:
    """Vertexwise closest-point thickness fields on CRUISE surfaces."""

    inner_surface: SurfaceMesh
    central_surface: SurfaceMesh | None
    mean_mm: float
    standard_deviation_mm: float
    nonzero_vertex_count: int


def levelset_subtraction_thickness(
    inner: nib.spatialimages.SpatialImage,
    outer: nib.spatialimages.SpatialImage,
    *,
    iso_level: float = 0.0,
) -> ThicknessVolumeResult:
    """Port the default ``GetThicknessCRUISE`` level-set subtraction mode.

    Voxel outputs are ``abs(inner - outer)`` on the inner grid. The historical
    summary statistics include only voxels whose two level sets lie on
    opposite sides of ``iso_level`` and use MATLAB's sample standard deviation.
    """
    inner_ras = _validated_thickness_image(inner, "inner level set")
    outer_ras = _validated_thickness_image(outer, "outer level set")
    if inner_ras.shape != outer_ras.shape or not np.allclose(
        inner_ras.affine, outer_ras.affine, rtol=0, atol=1e-4
    ):
        raise ValueError("Inner and outer level sets must share the same RAS grid")
    if not np.isfinite(iso_level):
        raise ValueError("Thickness iso-level must be finite")

    inner_values = np.asarray(inner_ras.dataobj, dtype=np.float32)
    outer_values = np.asarray(outer_ras.dataobj, dtype=np.float32)
    if not np.isfinite(inner_values).all() or not np.isfinite(outer_values).all():
        raise ValueError("Thickness level sets must contain only finite values")
    thickness = np.abs(inner_values - outer_values).astype(np.float32, copy=False)
    crosses = ((outer_values < iso_level) & (inner_values >= iso_level)) | (
        (outer_values >= iso_level) & (inner_values < iso_level)
    )
    sample = thickness[crosses].astype(np.float64)
    count = int(sample.size)
    mean = float(sample.mean()) if count else float("nan")
    stdev = float(sample.std(ddof=1)) if count > 1 else (0.0 if count == 1 else float("nan"))
    minimum = float(sample.min()) if count else float("nan")
    maximum = float(sample.max()) if count else float("nan")

    header = inner_ras.header.copy()
    header.set_data_dtype(np.float32)
    header.set_slope_inter(None, None)
    image = nib.Nifti1Image(thickness, inner_ras.affine, header=header)
    qform, qcode = inner_ras.get_qform(coded=True)
    _, scode = inner_ras.get_sform(coded=True)
    image.set_qform(qform if qcode else None, code=int(qcode), strip_shears=True)
    image.set_sform(inner_ras.affine, code=int(scode) if scode else 1)
    return ThicknessVolumeResult(image, mean, stdev, minimum, maximum, count)


def _validated_thickness_image(
    image: nib.spatialimages.SpatialImage,
    name: str,
) -> nib.spatialimages.SpatialImage:
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError(f"{name} must be a NIfTI image")
    if len(image.shape) != 3:
        raise ValueError(f"{name} must be 3D")
    return normalize_image_to_ras(image, source_name=name)


def _closest_points_on_triangles(
    point: np.ndarray,
    triangles: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return exact closest point and squared distance to candidate triangles."""
    a, b, c = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    ab = b - a
    ac = c - a
    normal = np.cross(ab, ac)
    normal_squared = np.einsum("ij,ij->i", normal, normal)
    projected = point - (
        np.einsum("ij,ij->i", point - a, normal)
        / np.maximum(normal_squared, np.finfo(float).tiny)
    )[:, None] * normal
    ap = projected - a
    d00 = np.einsum("ij,ij->i", ab, ab)
    d01 = np.einsum("ij,ij->i", ab, ac)
    d11 = np.einsum("ij,ij->i", ac, ac)
    d20 = np.einsum("ij,ij->i", ap, ab)
    d21 = np.einsum("ij,ij->i", ap, ac)
    denominator = d00 * d11 - d01 * d01
    v = np.divide(d11 * d20 - d01 * d21, denominator, out=np.zeros_like(d20), where=denominator > 0)
    w = np.divide(d00 * d21 - d01 * d20, denominator, out=np.zeros_like(d21), where=denominator > 0)
    inside = (v >= 0) & (w >= 0) & (v + w <= 1)
    candidates: list[np.ndarray] = []
    distances: list[np.ndarray] = []
    if np.any(inside):
        candidates.append(projected)
        distances.append(np.where(inside, np.einsum("ij,ij->i", projected - point, projected - point), np.inf))
    for start, end in ((a, b), (b, c), (c, a)):
        edge = end - start
        t = np.einsum("ij,ij->i", point - start, edge)
        t = np.divide(t, np.einsum("ij,ij->i", edge, edge), out=np.zeros_like(t), where=np.einsum("ij,ij->i", edge, edge) > 0)
        closest = start + np.clip(t, 0, 1)[:, None] * edge
        candidates.append(closest)
        distances.append(np.einsum("ij,ij->i", closest - point, closest - point))
    all_distances = np.stack(distances)
    rows = np.argmin(all_distances, axis=0)
    columns = np.arange(len(triangles))
    all_points = np.stack(candidates)
    return all_points[rows, columns], all_distances[rows, columns]


def _nearest_surface_points(
    points: np.ndarray,
    mesh: SurfaceMesh,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run the source-faithful JIST triangle tree on RAS millimetre points."""
    return JistSurfaceTree(mesh.vertices, mesh.faces).query_many(points)


def closest_point_thickness(
    inner: SurfaceMesh,
    outer: SurfaceMesh,
    *,
    central: SurfaceMesh | None = None,
) -> SurfaceThicknessResult:
    """Port ``ClosestPointThicknessForInner`` surface measurements.

    Inner point data contain inner-to-outer triangle distance, return distance
    from that closest outer point back to the inner surface, and their mean.
    If a central surface is supplied, its point data also project the mean
    thickness from the nearest inner triangle vertex, matching the source.
    """
    inner_tree = JistSurfaceTree(inner.vertices, inner.faces)
    outer_tree = JistSurfaceTree(outer.vertices, outer.faces)
    inner_to_outer_points, distance_inner_to_outer, _ = outer_tree.query_many(
        inner.vertices
    )
    _return_points, distance_outer_to_inner, _ = inner_tree.query_many(
        inner_to_outer_points
    )
    mean_thickness = (distance_inner_to_outer + distance_outer_to_inner) / 2.0
    inner_data = dict(inner.point_data)
    inner_data.update(
        {
            "distance_inner_to_outer_mm": distance_inner_to_outer,
            "distance_outer_to_inner_mm": distance_outer_to_inner,
            "mean_thickness_mm": mean_thickness,
        }
    )
    inner_with_data = SurfaceMesh(inner.vertices, inner.faces, point_data=inner_data)

    central_with_data: SurfaceMesh | None = None
    if central is not None:
        closest_on_inner, distance_to_inner, triangle_ids = inner_tree.query_many(
            central.vertices
        )
        triangle_vertices = inner.faces[triangle_ids]
        distances_to_vertices = np.linalg.norm(
            central.vertices[:, None, :] - inner.vertices[triangle_vertices], axis=2
        )
        nearest_corners = np.argmin(distances_to_vertices, axis=1)
        nearest_inner_vertex = triangle_vertices[np.arange(len(triangle_ids)), nearest_corners]
        central_data = dict(central.point_data)
        central_data.update(
            {
                "distance_inner_to_outer_mm": mean_thickness[nearest_inner_vertex],
                "distance_central_to_inner_mm": distance_to_inner,
                "distance_to_inner_vertex_mm": distances_to_vertices[
                    np.arange(len(triangle_ids)), nearest_corners
                ],
                "nearest_inner_vertex_index": nearest_inner_vertex.astype(np.float64),
            }
        )
        central_with_data = SurfaceMesh(
            central.vertices, central.faces, point_data=central_data
        )

    nonzero = mean_thickness[mean_thickness != 0]
    count = int(nonzero.size)
    mean = float(nonzero.mean()) if count else float("nan")
    stdev = float(nonzero.std(ddof=1)) if count > 1 else (0.0 if count == 1 else float("nan"))
    return SurfaceThicknessResult(
        inner_with_data,
        central_with_data,
        mean,
        stdev,
        count,
    )
