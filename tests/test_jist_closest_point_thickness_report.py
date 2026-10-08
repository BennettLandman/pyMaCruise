import json
from pathlib import Path


def test_selected_closest_point_thickness_report_keeps_tie_deviation_visible():
    report_path = (
        Path(__file__).resolve().parents[1]
        / "docs/reports/jist-closest-point-thickness-java-comparison-2026-10-07.json"
    )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    metrics = report["metrics"]
    assert report["synthetic_only"] is True
    assert metrics["comparison_pass"] is False
    assert metrics["max_absolute_field_difference"]["inner.distance_outer_to_inner_mm"] > 0.001
    assert metrics["max_absolute_field_difference"]["central.distance_inner_to_outer_mm"] > 0.001
    assert metrics["central_nearest_inner_vertex_index_mismatches"] == 21
    tie = report["tie_diagnostic"]
    assert tie["java_selected_outer_face"] != tie["python_selected_outer_face"]
    assert abs(tie["java_outer_distance_mm"] - tie["python_outer_distance_mm"]) < 1e-6
    assert abs(
        tie["java_return_distance_mm"] - tie["python_return_distance_from_java_hit_mm"]
    ) < 1e-6
