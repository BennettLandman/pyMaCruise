"""Pin the full-grid selected-source JIST surface comparison evidence."""
from __future__ import annotations

import json
from pathlib import Path


REPORT = (
    Path(__file__).parent.parent
    / "docs/reports/current-base-jist-surface-comparison-2026-10-07.json"
)


def test_brain_resolution_jist_surface_report_passes_p08_and_retains_topology_gap() -> None:
    report = json.loads(REPORT.read_text())
    assert report["synthetic_only"] is True
    assert report["evidence_type"] == "synthetic selected-source comparison"
    assert report["source"]["method"].endswith("18, 0.0f, false)")
    assert report["fixture"]["shape"] == [172, 216, 164]
    assert report["fixture"]["orientation"] == "RAS"

    for surface in ("inner", "central", "outer"):
        metrics = report["metrics"]["per_surface"][surface]
        position = metrics["vertex_to_triangle_surface_distance_mm"]
        assert position["symmetric_mean"] <= 0.01
        assert position["symmetric_max"] <= 1.0
        assert metrics["jist_topology"]["connected_components"] == 1
        assert metrics["lewiner_topology"]["connected_components"] > 1
        assert (
            metrics["triangle_connectivity_after_nearest_vertex_mapping"][
                "matching_fraction_of_jist_faces"
            ]
            < 0.70
        )

        lorensen = report["metrics"]["lorensen_per_surface"][surface]
        lorensen_position = lorensen["vertex_to_triangle_surface_distance_mm"]
        assert lorensen_position["symmetric_mean"] <= 0.01
        assert lorensen_position["symmetric_max"] <= 1.0
        assert lorensen["jist_topology"]["connected_components"] == 1
        assert lorensen["lewiner_topology"]["connected_components"] > 1
        assert (
            lorensen["triangle_connectivity_after_nearest_vertex_mapping"][
                "matching_fraction_of_jist_faces"
            ]
            < 0.70
        )

    assert [
        report["metrics"]["lorensen_per_surface"][name]["lewiner_topology"][
            "connected_components"
        ]
        for name in ("inner", "central", "outer")
    ] == [2, 2, 28]

    thickness = report["metrics"]["thickness"]
    assert abs(thickness["mean_delta_mm"] - (-0.00008807700348656411)) < 1e-12
    assert abs(thickness["lorensen_mean_delta_mm"] - 0.0000124834290862) < 1e-12
