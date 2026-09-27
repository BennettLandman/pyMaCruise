"""Synthetic-first comparators for image, segmentation, and surface outputs."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

def _same_shape(reference: np.ndarray, candidate: np.ndarray) -> None:
    if reference.shape != candidate.shape:
        raise ValueError(
            f"Array shapes differ: {reference.shape} versus {candidate.shape}"
        )


def dice_coefficient(reference: np.ndarray, candidate: np.ndarray) -> float:
    """Return binary Dice; both-empty inputs score 1 and one-empty scores 0."""
    _same_shape(reference, candidate)
    ref = np.asarray(reference, dtype=bool)
    cand = np.asarray(candidate, dtype=bool)
    denominator = int(ref.sum()) + int(cand.sum())
    if denominator == 0:
        return 1.0
    return 2.0 * float(np.count_nonzero(ref & cand)) / denominator


def jaccard_coefficient(reference: np.ndarray, candidate: np.ndarray) -> float:
    """Return binary intersection-over-union; both-empty inputs score 1."""
    _same_shape(reference, candidate)
    ref = np.asarray(reference, dtype=bool)
    cand = np.asarray(candidate, dtype=bool)
    union = int(np.count_nonzero(ref | cand))
    if union == 0:
        return 1.0
    return float(np.count_nonzero(ref & cand)) / union


def exact_voxel_metrics(reference: np.ndarray, candidate: np.ndarray) -> dict[str, int | bool]:
    """Return exact-equality status and differing voxel count."""
    _same_shape(reference, candidate)
    differences = int(np.count_nonzero(np.asarray(reference) != np.asarray(candidate)))
    return {
        "exact": differences == 0,
        "n_different": differences,
        "n_values": int(reference.size),
    }


def label_dice_scores(
    reference: np.ndarray,
    candidate: np.ndarray,
    *,
    labels: tuple[int, ...] | None = None,
    background: int = 0,
) -> dict[int, float]:
    """Return per-label Dice scores, omitting background by default."""
    _same_shape(reference, candidate)
    if labels is None:
        labels = tuple(
            int(label)
            for label in np.union1d(np.unique(reference), np.unique(candidate))
            if label != background
        )
    return {
        int(label): dice_coefficient(reference == label, candidate == label)
        for label in labels
    }


def voxel_error_metrics(reference: np.ndarray, candidate: np.ndarray) -> dict[str, float | int]:
    """Return MAE, RMSE, and maximum absolute error for equal-shape arrays."""
    _same_shape(reference, candidate)
    ref = np.asarray(reference, dtype=np.float64)
    cand = np.asarray(candidate, dtype=np.float64)
    if not np.isfinite(ref).all() or not np.isfinite(cand).all():
        raise ValueError("Voxel comparison inputs must contain only finite values")
    delta = np.abs(ref - cand)
    return {
        "mae": float(delta.mean()),
        "rmse": float(np.sqrt(np.mean(np.square(delta)))),
        "max_abs": float(delta.max(initial=0.0)),
        "n_values": int(delta.size),
    }


def point_cloud_distance_metrics(
    reference_ras_mm: np.ndarray,
    candidate_ras_mm: np.ndarray,
) -> dict[str, float | int]:
    """Compare finite 3D RAS+ point samples with symmetric nearest neighbors.

    Distances are Euclidean millimetres. This is a point-sample comparator;
    triangle-surface distances require sufficiently dense, documented samples
    or a geometry-aware closest-point implementation.
    """
    reference = np.asarray(reference_ras_mm, dtype=np.float64)
    candidate = np.asarray(candidate_ras_mm, dtype=np.float64)
    for name, points in (("reference", reference), ("candidate", candidate)):
        if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] == 0:
            raise ValueError(f"{name} points must be a non-empty N×3 array")
        if not np.isfinite(points).all():
            raise ValueError(f"{name} points must contain only finite RAS coordinates")

    ref_to_cand = cKDTree(candidate).query(reference, k=1)[0]
    cand_to_ref = cKDTree(reference).query(candidate, k=1)[0]
    return {
        "hausdorff_mm": float(max(ref_to_cand.max(), cand_to_ref.max())),
        "mean_absolute_distance_mm": float(
            (ref_to_cand.mean() + cand_to_ref.mean()) / 2.0
        ),
        "reference_to_candidate_mean_mm": float(ref_to_cand.mean()),
        "candidate_to_reference_mean_mm": float(cand_to_ref.mean()),
        "n_reference_points": int(reference.shape[0]),
        "n_candidate_points": int(candidate.shape[0]),
    }


def transform_error_metrics(
    reference_ras: np.ndarray,
    candidate_ras: np.ndarray,
) -> dict[str, float]:
    """Compare two homogeneous RAS transforms using rigid and matrix errors."""
    reference = np.asarray(reference_ras, dtype=np.float64)
    candidate = np.asarray(candidate_ras, dtype=np.float64)
    for name, transform in (("reference", reference), ("candidate", candidate)):
        if transform.shape != (4, 4) or not np.isfinite(transform).all():
            raise ValueError(f"{name} transform must be a finite 4×4 matrix")
        if not np.allclose(transform[3], (0.0, 0.0, 0.0, 1.0), atol=1e-8):
            raise ValueError(f"{name} transform must have homogeneous final row")
        rotation = transform[:3, :3]
        if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-6) or not np.isclose(
            np.linalg.det(rotation), 1.0, atol=1e-6
        ):
            raise ValueError(f"{name} transform must contain a proper rigid rotation")
    relative_rotation = reference[:3, :3].T @ candidate[:3, :3]
    cosine = (np.trace(relative_rotation) - 1.0) / 2.0
    angle = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    translation_mm = np.linalg.norm(reference[:3, 3] - candidate[:3, 3])
    return {
        "translation_error_mm": float(translation_mm),
        "rotation_error_degrees": float(angle),
        "max_matrix_abs_error": float(np.max(np.abs(reference - candidate))),
    }


def _point_to_triangles_distance(
    point: np.ndarray, triangles: np.ndarray
) -> np.ndarray:
    """Compute exact Euclidean distance from one point to each triangle."""
    a, b, c = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    ab = b - a
    ac = c - a
    normal = np.cross(ab, ac)
    normal_sq = np.einsum("ij,ij->i", normal, normal)
    ap = point - a

    # Orthogonal projection onto each triangle plane; retain it only when its
    # barycentric coordinates place it inside the triangle.
    plane_offset = np.divide(
        np.einsum("ij,ij->i", ap, normal),
        normal_sq,
        out=np.zeros_like(normal_sq),
        where=normal_sq > 1e-20,
    )
    projected = point - plane_offset[:, None] * normal
    v2 = projected - a
    d00 = np.einsum("ij,ij->i", ab, ab)
    d01 = np.einsum("ij,ij->i", ab, ac)
    d11 = np.einsum("ij,ij->i", ac, ac)
    d20 = np.einsum("ij,ij->i", v2, ab)
    d21 = np.einsum("ij,ij->i", v2, ac)
    denominator = d00 * d11 - d01 * d01
    v = np.divide(
        d11 * d20 - d01 * d21,
        denominator,
        out=np.full_like(denominator, -1.0),
        where=denominator > 1e-20,
    )
    w = np.divide(
        d00 * d21 - d01 * d20,
        denominator,
        out=np.full_like(denominator, -1.0),
        where=denominator > 1e-20,
    )
    inside = (v >= 0) & (w >= 0) & (v + w <= 1) & (normal_sq > 1e-20)
    distances_sq = np.full(triangles.shape[0], np.inf, dtype=np.float64)
    distances_sq[inside] = np.square(plane_offset[inside]) * normal_sq[inside]

    # A closest point outside the face lies on one of its three edges. This
    # also handles degenerate triangles through their non-zero edges.
    for start, end in ((a, b), (b, c), (c, a)):
        edge = end - start
        edge_sq = np.einsum("ij,ij->i", edge, edge)
        fraction = np.divide(
            np.einsum("ij,ij->i", point - start, edge),
            edge_sq,
            out=np.zeros_like(edge_sq),
            where=edge_sq > 1e-20,
        )
        fraction = np.clip(fraction, 0.0, 1.0)
        closest = start + fraction[:, None] * edge
        edge_delta = point - closest
        edge_distance_sq = np.einsum("ij,ij->i", edge_delta, edge_delta)
        distances_sq = np.minimum(distances_sq, edge_distance_sq)
    return np.sqrt(distances_sq)


def _vertices_to_surface_distances(
    query_vertices: np.ndarray,
    target_vertices: np.ndarray,
    target_faces: np.ndarray,
) -> np.ndarray:
    triangles = target_vertices[target_faces]
    centroids = triangles.mean(axis=1)
    radii = np.linalg.norm(triangles - centroids[:, None, :], axis=2).max(axis=1)
    max_radius = float(radii.max())
    tree = cKDTree(centroids)
    distances = np.empty(query_vertices.shape[0], dtype=np.float64)
    for index, point in enumerate(query_vertices):
        initial_distance, _ = tree.query(point, k=1)
        candidate_ids = tree.query_ball_point(
            point, float(initial_distance) + max_radius + 1e-10
        )
        candidate_triangles = triangles[np.asarray(candidate_ids, dtype=np.int64)]
        distances[index] = _point_to_triangles_distance(point, candidate_triangles).min()
    return distances


def mesh_surface_distance_metrics(
    reference_vertices_ras_mm: np.ndarray,
    reference_faces: np.ndarray,
    candidate_vertices_ras_mm: np.ndarray,
    candidate_faces: np.ndarray,
) -> dict[str, float | int]:
    """Measure symmetric mesh-vertex to continuous-triangle-surface distances.

    Closest points are computed on triangle interiors and edges, not just on
    the other mesh's vertices. Hausdorff and mean distances use each mesh's
    vertices as deterministic surface samples; units are RAS+ millimetres.
    """
    reference_vertices = np.asarray(reference_vertices_ras_mm, dtype=np.float64)
    candidate_vertices = np.asarray(candidate_vertices_ras_mm, dtype=np.float64)
    reference_triangles = np.asarray(reference_faces, dtype=np.int64)
    candidate_triangles = np.asarray(candidate_faces, dtype=np.int64)
    for name, vertices, faces in (
        ("reference", reference_vertices, reference_triangles),
        ("candidate", candidate_vertices, candidate_triangles),
    ):
        if vertices.ndim != 2 or vertices.shape[1] != 3 or vertices.shape[0] == 0:
            raise ValueError(f"{name} vertices must be a non-empty N×3 array")
        if faces.ndim != 2 or faces.shape[1] != 3 or faces.shape[0] == 0:
            raise ValueError(f"{name} faces must be a non-empty M×3 array")
        if not np.isfinite(vertices).all():
            raise ValueError(f"{name} vertices must contain only finite RAS coordinates")
        if faces.min() < 0 or faces.max() >= vertices.shape[0]:
            raise ValueError(f"{name} face indices are outside the vertex array")

    ref_to_cand = _vertices_to_surface_distances(
        reference_vertices, candidate_vertices, candidate_triangles
    )
    cand_to_ref = _vertices_to_surface_distances(
        candidate_vertices, reference_vertices, reference_triangles
    )
    return {
        "hausdorff_mm": float(max(ref_to_cand.max(), cand_to_ref.max())),
        "mean_absolute_distance_mm": float((ref_to_cand.mean() + cand_to_ref.mean()) / 2),
        "reference_to_candidate_mean_mm": float(ref_to_cand.mean()),
        "candidate_to_reference_mean_mm": float(cand_to_ref.mean()),
        "n_reference_vertices": int(reference_vertices.shape[0]),
        "n_candidate_vertices": int(candidate_vertices.shape[0]),
    }


def write_metrics_report(path: str | Path, metrics: dict[str, Any]) -> None:
    """Write nested comparison metrics as strict, deterministic JSON."""
    destination = Path(path)
    if destination.suffix.lower() != ".json":
        raise ValueError("Validation metrics reports must use the .json suffix")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(metrics, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def mesh_topology_metrics(
    vertices_ras_mm: np.ndarray, faces: np.ndarray
) -> dict[str, int | bool]:
    """Return basic triangle-mesh validity and connected-component metrics."""
    vertices = np.asarray(vertices_ras_mm, dtype=np.float64)
    triangles = np.asarray(faces, dtype=np.int64)
    if vertices.ndim != 2 or vertices.shape[1] != 3 or vertices.shape[0] == 0:
        raise ValueError("Mesh vertices must be a non-empty N×3 array")
    if triangles.ndim != 2 or triangles.shape[1] != 3 or triangles.shape[0] == 0:
        raise ValueError("Mesh faces must be a non-empty M×3 array")
    if not np.isfinite(vertices).all():
        raise ValueError("Mesh vertices must contain only finite RAS coordinates")
    if triangles.min() < 0 or triangles.max() >= vertices.shape[0]:
        raise ValueError("Mesh face indices are outside the vertex array")

    edges = np.sort(
        np.concatenate(
            (triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]])
        ),
        axis=1,
    )
    unique_edges, edge_counts = np.unique(edges, axis=0, return_counts=True)
    edge_histogram = Counter(edge_counts.tolist())
    edge_start, edge_end = unique_edges.T
    adjacency_rows = np.concatenate((edge_start, edge_end))
    adjacency_cols = np.concatenate((edge_end, edge_start))
    adjacency = sparse.coo_matrix(
        (
            np.ones(adjacency_rows.size, dtype=np.uint8),
            (adjacency_rows, adjacency_cols),
        ),
        shape=(vertices.shape[0], vertices.shape[0]),
    ).tocsr()
    component_count, _ = connected_components(adjacency, directed=False)
    return {
        "n_vertices": int(vertices.shape[0]),
        "n_faces": int(triangles.shape[0]),
        "n_edges": int(unique_edges.shape[0]),
        "n_components": int(component_count),
        "n_open_edges": int(edge_histogram.get(1, 0)),
        "n_nonmanifold_edges": int(
            sum(count for uses, count in edge_histogram.items() if uses > 2)
        ),
        "watertight_by_edge_incidence": bool(np.all(edge_counts == 2)),
    }
