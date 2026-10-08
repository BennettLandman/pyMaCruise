"""Exercise every packaged cube case captured from the selected JIST source."""
from __future__ import annotations

import json
from importlib.resources import files

import numpy as np
import pytest

from pymacruise.surfaces import _jist_connectivity18_cases, _jist_connectivity18_surface


TABLE_PATH = files("pymacruise").joinpath(
    "data/jist_isosurface_connectivity18_cases.json"
)


def test_packaged_table_records_selected_java_source_and_full_case_coverage() -> None:
    table = json.loads(TABLE_PATH.read_text(encoding="utf-8"))
    assert table["source_method"] == "solveOriginal(ImageDataFloat,18,0,false)"
    assert table["source_revision"] == "CVS r1.7"
    assert len(table["patterns"]) == 256
    assert table["validation"].startswith("all 256 sign patterns; five distinct magnitude assignments")
    assert len(_jist_connectivity18_cases()) == 256


@pytest.mark.parametrize("pattern", range(256))
def test_cube_output_matches_java_derived_case_template(pattern: int) -> None:
    """Check runtime interpolation and indexed faces for each Java case."""
    data = np.empty((2, 2, 2), dtype=np.float32)
    corner_values = np.empty(8, dtype=np.float32)
    coords = np.asarray(
        [(bit & 1, (bit >> 1) & 1, (bit >> 2) & 1) for bit in range(8)],
        dtype=np.float32,
    )
    for bit, (x, y, z) in enumerate(coords.astype(int)):
        magnitude = np.float32(1.0 + bit * 0.137)
        value = magnitude if (pattern >> bit) & 1 else -magnitude
        data[x, y, z] = value
        corner_values[bit] = value

    if pattern in (0, 255):
        with pytest.raises(ValueError, match="does not cross"):
            _jist_connectivity18_surface(data)
        return

    actual_vertices, actual_faces = _jist_connectivity18_surface(data)
    case = _jist_connectivity18_cases()[pattern]
    edge_pairs = tuple(
        (a, b)
        for a in range(8)
        for b in range(a + 1, 8)
        if ((a ^ b) & ((a ^ b) - 1)) == 0
    )
    expected_vertices = []
    for edge_id in case["edges"]:
        a, b = edge_pairs[edge_id]
        fraction = -float(corner_values[a]) / float(corner_values[b] - corner_values[a])
        expected_vertices.append(
            (coords[a] + fraction * (coords[b] - coords[a])).astype(np.float32)
        )
    np.testing.assert_array_equal(actual_faces, np.asarray(case["faces"], dtype=np.int64))
    np.testing.assert_allclose(
        actual_vertices,
        np.asarray(expected_vertices, dtype=np.float32),
        rtol=0,
        atol=1e-7,
    )
