"""Pin the selected-source JIST point-to-triangle comparison evidence."""
from __future__ import annotations

import json
from pathlib import Path


REPORT = (
    Path(__file__).parent.parent
    / "docs/reports/jist-triangle-distance-java-comparison-2026-10-07.json"
)


def test_jist_triangle_distance_matches_on_synthetic_queries() -> None:
    report = json.loads(REPORT.read_text())
    assert report["synthetic_only"] is True
    assert report["fixture"]["query_count"] == 256
    assert report["fixture"]["triangle_count"] > 1000
    assert report["metrics"]["all_queries_pass"] is True
    assert report["metrics"]["distance_absolute_error_mm"]["max"] < 1e-5
    assert report["metrics"]["closest_point_euclidean_error_mm"]["max"] < 1e-3
    assert "does not execute the selected SurfaceIntersector KD traversal" in report["method_scope"]
