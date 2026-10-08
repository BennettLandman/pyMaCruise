import json
from pathlib import Path


def test_selected_closest_point_thickness_report_records_java_match_and_baseline_ties():
    report_path = (
        Path(__file__).resolve().parents[1]
        / "docs/reports/jist-closest-point-thickness-java-comparison-2026-10-07.json"
    )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    metrics = report["metrics"]
    assert report["synthetic_only"] is True
    assert metrics["comparison_pass"] is True
    assert metrics["max_absolute_field_difference"][
        "central.distance_to_inner_vertex_mm"
    ] < 1e-6
    assert metrics["central_nearest_inner_vertex_index_mismatches"] == 0
    ties = report["centroid_tree_baseline"]["metrics"]
    assert ties["inner_to_outer_face_mismatches"] == 569
    assert ties["central_to_inner_face_mismatches"] == 932
    assert ties["max_outer_to_inner_distance_error_mm"] > 0.01
    faces = report["face_order_diagnostics"]
    assert faces["inner"]["query_count"] == 2194
    assert faces["inner"]["java_python_selected_face_mismatches"] == 569
    assert faces["central"]["query_count"] == 3062
    assert faces["central"]["java_python_selected_face_mismatches"] == 932
    assert faces["inner"]["mismatches_where_java_face_is_highest_face_within_1e-6_mm_tie"] < (
        faces["inner"]["tie_policy_sample_size"] * 0.6
    )
    port = report["python_jist_tree_port"]
    assert port["status"].startswith("connected to production")
    assert all(tree["tree_structure_exact"] for tree in port["tree_builds"].values())
    assert port["metrics"]["central_nearest_vertex_index_mismatches"] == 0
    assert port["metrics"]["max_outer_to_inner_distance_error_mm"] == 0


def test_brain_resolution_selected_thickness_wrapper_matches_java_exactly():
    report_path = (
        Path(__file__).resolve().parents[1]
        / "docs/reports/jist-closest-point-thickness-brain-resolution-java-comparison-2026-10-07.json"
    )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["synthetic_only"] is True
    assert report["metrics"]["comparison_pass"] is True
    assert report["metrics"]["max_absolute_field_difference"] == {
        "central.distance_central_to_inner_mm": 0.0,
        "central.distance_inner_to_outer_mm": 0.0,
        "central.distance_to_inner_vertex_mm": 0.0,
        "central.nearest_inner_vertex_index": 0.0,
        "inner.distance_inner_to_outer_mm": 0.0,
        "inner.distance_outer_to_inner_mm": 0.0,
        "inner.mean_thickness_mm": 0.0,
    }
    assert report["metrics"]["central_nearest_inner_vertex_index_mismatches"] == 0
    assert report["fixture"]["mesh_counts"]["inner"]["faces"] == 232328
    assert report["fixture"]["mesh_counts"]["central"]["vertices"] == 121074
    assert report["runtime"]["python_process_peak_rss_gb"] < 2.0
