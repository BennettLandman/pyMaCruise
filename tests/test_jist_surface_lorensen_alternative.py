from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).parent
SCRIPT = ROOT.parent / "scripts" / "compare_jist_surface_lorensen.py"


def test_lorensen_alternative_meets_sampled_surface_gates_but_keeps_connectivity_gap() -> None:
    result = subprocess.run(
        [sys.executable, str(SCRIPT)],
        check=True,
        capture_output=True,
        text=True,
        cwd=ROOT.parent,
    )
    comparison = json.loads(result.stdout)

    assert comparison["method"] == "scikit-image marching_cubes(method='lorensen')"
    assert set(comparison["per_surface"]) == {"inner", "central", "outer"}
    for surface in comparison["per_surface"].values():
        connectivity = surface["triangle_connectivity_after_nearest_vertex_mapping"]
        sampled_distance = surface["vertex_to_triangle_surface_distance_mm"]
        assert connectivity["vertex_mapping_is_bijective"]
        assert connectivity["jist_face_count"] == surface["jist_topology"]["faces"]
        assert connectivity["lewiner_face_count"] == surface["lewiner_topology"]["faces"]
        assert connectivity["matching_fraction_of_jist_faces"] < 0.8
        # P08 defines symmetric mean nearest-surface distance from vertices
        # projected to the opposite mesh's triangles, with a 0.01 mm limit.
        assert sampled_distance["symmetric_mean"] <= 0.01
        # This bounds the same deterministic sample set; it is not a
        # continuous triangle-to-triangle Hausdorff certificate.
        assert sampled_distance["symmetric_max"] <= 1.0
        for topology in (surface["jist_topology"], surface["lewiner_topology"]):
            assert topology["boundary_edges"] == 0
            assert topology["nonmanifold_edges"] == 0
            assert topology["connected_components"] == 1
            assert topology["euler_characteristic"] == 2
