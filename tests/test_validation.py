from __future__ import annotations

import numpy as np
import pytest
import json

from pymacruise.synthetic import _make_icosphere
from pymacruise.validation import (
    dice_coefficient,
    exact_voxel_metrics,
    jaccard_coefficient,
    mesh_surface_distance_metrics,
    label_dice_scores,
    mesh_topology_metrics,
    point_cloud_distance_metrics,
    transform_error_metrics,
    voxel_error_metrics,
    write_metrics_report,
)


def test_dice_and_label_dice_edge_cases() -> None:
    assert dice_coefficient(np.zeros(4), np.zeros(4)) == 1.0
    assert dice_coefficient(np.ones(4), np.zeros(4)) == 0.0
    assert jaccard_coefficient(np.zeros(4), np.zeros(4)) == 1.0
    reference = np.array([0, 1, 1, 2, 2, 0])
    candidate = np.array([0, 1, 2, 2, 0, 0])
    scores = label_dice_scores(reference, candidate)
    assert scores == {1: pytest.approx(2 / 3), 2: pytest.approx(0.5)}
    assert jaccard_coefficient(reference == 1, candidate == 1) == pytest.approx(0.5)
    assert exact_voxel_metrics(reference, reference.copy()) == {
        "exact": True,
        "n_different": 0,
        "n_values": 6,
    }


def test_voxel_error_metrics_and_perturbed_output() -> None:
    reference = np.zeros((2, 2), dtype=np.float32)
    candidate = reference.copy()
    candidate[0, 0] = 0.5
    metrics = voxel_error_metrics(reference, candidate)
    assert metrics["mae"] == pytest.approx(0.125)
    assert metrics["rmse"] == pytest.approx(0.25)
    assert metrics["max_abs"] == pytest.approx(0.5)
    assert metrics["max_abs"] > 0.01


def test_symmetric_point_distances_report_millimetres() -> None:
    reference = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    candidate = reference + np.array([0.02, 0.0, 0.0])
    metrics = point_cloud_distance_metrics(reference, candidate)
    assert metrics["hausdorff_mm"] == pytest.approx(0.02)
    assert metrics["mean_absolute_distance_mm"] == pytest.approx(0.02)
    assert metrics["hausdorff_mm"] > 0.01


def test_mesh_surface_distance_handles_different_tessellations_and_offsets() -> None:
    square = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    coarse_faces = np.array([[0, 1, 2], [0, 2, 3]])
    fine = np.vstack((square, np.array([[0.5, 0.5, 0.0]])))
    fine_faces = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])

    exact = mesh_surface_distance_metrics(square, coarse_faces, fine, fine_faces)
    assert exact["hausdorff_mm"] == pytest.approx(0.0)
    assert exact["mean_absolute_distance_mm"] <= 0.01

    shifted = fine.copy()
    shifted[:, 2] = 0.02
    outside = mesh_surface_distance_metrics(square, coarse_faces, shifted, fine_faces)
    assert outside["hausdorff_mm"] == pytest.approx(0.02)
    assert outside["mean_absolute_distance_mm"] == pytest.approx(0.02)
    assert outside["mean_absolute_distance_mm"] > 0.01


def test_transform_errors_and_machine_readable_report(tmp_path) -> None:
    reference = np.eye(4)
    candidate = np.eye(4)
    candidate[:3, 3] = (0.0, 0.3, 0.4)
    metrics = transform_error_metrics(reference, candidate)
    assert metrics["translation_error_mm"] == pytest.approx(0.5)
    assert metrics["rotation_error_degrees"] == pytest.approx(0.0)

    report = tmp_path / "metrics" / "report.json"
    write_metrics_report(report, {"transform": metrics, "passed": False})
    assert json.loads(report.read_text()) == {"passed": False, "transform": metrics}
    with pytest.raises(ValueError, match=".json"):
        write_metrics_report(tmp_path / "report.txt", metrics)


def test_mesh_topology_metrics_detect_closed_and_perturbed_meshes() -> None:
    vertices, faces = _make_icosphere(subdivisions=1)
    metrics = mesh_topology_metrics(vertices, faces)
    assert metrics["watertight_by_edge_incidence"] is True
    assert metrics["n_components"] == 1
    assert metrics["n_open_edges"] == 0

    broken = faces[:-1]
    broken_metrics = mesh_topology_metrics(vertices, broken)
    assert broken_metrics["watertight_by_edge_incidence"] is False
    assert broken_metrics["n_open_edges"] > 0
