"""Pin selected-source SurfaceIntersector comparison evidence."""
from __future__ import annotations

import json
from pathlib import Path


REPORT = (
    Path(__file__).parent.parent
    / "docs/reports/jist-surface-intersector-java-comparison-2026-10-07.json"
)


def test_jist_surface_intersector_matches_python_on_synthetic_queries() -> None:
    report = json.loads(REPORT.read_text())
    assert report["synthetic_only"] is True
    assert report["fixture"]["query_count"] == 256
    assert report["metrics"]["all_queries_pass"] is True
    assert report["metrics"]["distance_absolute_error_mm"]["max"] < 1e-5
    assert report["metrics"]["closest_point_euclidean_error_mm"]["max"] < 1e-3
    assert "full ClosestPointThicknessForInner wrapper" in report["method_scope"]
    assert all(len(value) == 64 for value in report["sources"].values())
