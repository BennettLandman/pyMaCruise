# Selected JIST closest-point thickness on brain-resolution synthetic meshes

The selected VandyCRUISE thickness wrapper and production Python API ran on the same committed Java `IsoSurfaceOnGrid` meshes. This isolates the wrapper at the `(172,216,164)` base-core scale; upstream processing and mesh extraction are held fixed.

| Mesh | Vertices | Faces |
|---|---:|---:|
| inner | 116,166 | 232,328 |
| central | 121,074 | 242,144 |
| outer | 123,828 | 247,652 |

| Field | Maximum absolute difference | Mean absolute difference |
|---|---:|---:|
| `inner.distance_inner_to_outer_mm` | 0 | 0 |
| `inner.distance_outer_to_inner_mm` | 0 | 0 |
| `inner.mean_thickness_mm` | 0 | 0 |
| `central.distance_inner_to_outer_mm` | 0 | 0 |
| `central.distance_central_to_inner_mm` | 0 | 0 |
| `central.distance_to_inner_vertex_mm` | 0 | 0 |
| `central.nearest_inner_vertex_index` | 0 | 0 |

**Selected wrapper comparison: PASS.** Maximum physical-field error is `0 mm` against the `0.001 mm` scalar gate; central nearest-vertex index mismatches: 0.

The Python API took `473.8 s` and peaked at `0.625 GB` process RSS. Java ran with a 2 GiB heap cap. This is selected-source synthetic evidence, not whole-workflow equivalence or real-MRI validation. Mesh extraction connectivity differences remain open. See the [JSON report](jist-closest-point-thickness-brain-resolution-java-comparison-2026-10-07.json) for source and fixture hashes. Reproduce with `./.venv/bin/python scripts/compare_jist_surface_current_base.py` followed by `./.venv/bin/python scripts/compare_jist_closest_point_thickness_brain_resolution.py`.
