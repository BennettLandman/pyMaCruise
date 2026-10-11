# JIST `IsoSurfaceOnGrid` versus Lewiner: folded synthetic phantom

This report compares the selected Java `IsoSurfaceOnGrid.solveOriginal` implementation with PyMaCRUISE's scikit-image Lewiner extraction on the exact same three nested float32 level sets. The source wrapper selects connectivity 18, iso-level 0, and `inclusive=false`.

Source `IsoSurfaceOnGrid.java` SHA-256: `a323ac5eed8b66ad62a68cd93a6bc7bcf01281453fc6541ee5c4ccc529d40c38`. Wrapper SHA-256: `3f7b54e099f2550f42a572e96070642b8124a93ab1b45eca1aa544ecec1758b4`.

## Surface metrics

| Level set | JIST V/F/χ | Lewiner V/F/χ | JIST boundary/nonmanifold edges | Lewiner boundary/nonmanifold edges | Symmetric vertex NN mean / max (mm) | Symmetric vertex-to-triangle mean / max (mm) |
|---|---:|---:|---:|---:|---:|---:|
| inner | 2194/4384/2 | 2194/4384/2 | 0/0 | 0/0 | 0.000000 / 0.000001 | 0.000000 / 0.000000 |
| central | 3062/6120/2 | 3062/6120/2 | 0/0 | 0/0 | 0.000000 / 0.000002 | 0.000000 / 0.000000 |
| outer | 4070/8136/2 | 4070/8136/2 | 0/0 | 0/0 | 0.000000 / 0.000002 | 0.000000 / 0.000000 |

Vertex counts, topology invariants, and unordered triangle connectivity match on this phantom after bijectively mapping the near-identical Java vertices to Lewiner vertex indices. The vertex-to-triangle measurements are zero at reported precision; the vertex-nearest maximum remains a few micrometres because of tiny coordinate differences.

| Level set | Matching unordered triangles | Matching fraction |
|---|---:|---:|
| inner | 4384/4384 | 100.0% |
| central | 6120/6120 | 100.0% |
| outer | 8136/8136 | 100.0% |

## Downstream thickness impact

Using the same `closest_point_thickness` operator on the two nested mesh sets: JIST geometry mean thickness `3.558754 mm` (SD `0.367832`, 2194 nonzero inner vertices); Lewiner geometry mean thickness `3.558754 mm` (SD `0.367832`, 2194 nonzero inner vertices); difference `+0.000000 mm`.

## Scope and limitations

- The harness compiles and invokes the selected `IsoSurfaceOnGrid` source directly. Minimal API stubs replace JIST progress/image containers and normal generation; the harness records the Java algorithm's actual returned vertices and triangle indices.
- Both methods consume identical float32 values. Per-volume SHA-256 digests are checked by the runner. The fixture's identity affine makes index coordinates equal to RAS millimetres.
- Vertex-nearest distances depend on tessellation density. Vertex-to-triangle distances here sample mesh vertices against the other mesh's triangles; they are not a continuous triangle-to-triangle Hausdorff bound.
- All unordered faces match after vertex mapping for this fixture. This removes the observed synthetic triangulation-connectivity risk for this case, but does not establish matching behavior for other level sets or resolve the separate human ribbon-quality deficit.
- Thickness uses the same Python closest-point operator on both extracted mesh sets. This quantifies geometry's effect on that measure; it is not a direct Java thickness implementation comparison.
- This single synthetic folded fixture does not establish clinical or end-to-end equivalence. Matching this case would not remove the need for broader shapes, adversarial topology cases, and source-version verification.

Full metric output and source/input provenance: [`comparison.json`](../../tests/data/synthetic/jist_surface_comparison/comparison.json). Fixture generator: [`make_phantom.py`](../../tests/data/synthetic/jist_surface_comparison/make_phantom.py).
