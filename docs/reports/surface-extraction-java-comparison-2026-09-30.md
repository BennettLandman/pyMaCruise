# JIST `IsoSurfaceOnGrid` versus Lewiner: folded synthetic phantom

This report compares the selected Java `IsoSurfaceOnGrid.solveOriginal` implementation with PyMaCRUISE's scikit-image Lewiner extraction on the exact same three nested float32 level sets. The source wrapper selects connectivity 18, iso-level 0, and `inclusive=false`.

Source `IsoSurfaceOnGrid.java` SHA-256: `a323ac5eed8b66ad62a68cd93a6bc7bcf01281453fc6541ee5c4ccc529d40c38`. Wrapper SHA-256: `3f7b54e099f2550f42a572e96070642b8124a93ab1b45eca1aa544ecec1758b4`.

## Surface metrics

| Level set | JIST V/F/χ | Lewiner V/F/χ | JIST boundary/nonmanifold edges | Lewiner boundary/nonmanifold edges | Symmetric vertex NN mean / max (mm) | Symmetric vertex-to-triangle mean / max (mm) |
|---|---:|---:|---:|---:|---:|---:|
| inner | 2194/4384/2 | 2194/4384/2 | 0/0 | 0/0 | 0.000000 / 0.000001 | 0.000000 / 0.000001 |
| central | 3062/6120/2 | 3062/6120/2 | 0/0 | 0/0 | 0.000000 / 0.000002 | 0.000000 / 0.000002 |
| outer | 4070/8136/2 | 4070/8136/2 | 0/0 | 0/0 | 0.000000 / 0.000002 | 0.000000 / 0.000002 |

The vertex-count/topology invariants match on this phantom, but this does not mean the triangulations match. After bijectively mapping near-identical Java vertices to Lewiner vertex indices, unordered triangles match by the fractions below. The remaining face connectivity can change interpolated triangle interiors and mesh-dependent measurements.

| Level set | Matching unordered triangles | Matching fraction |
|---|---:|---:|
| inner | 3034/4384 | 69.2% |
| central | 4320/6120 | 70.6% |
| outer | 5640/8136 | 69.3% |

## Downstream thickness impact

Using the same `closest_point_thickness` operator on the two nested mesh sets: JIST geometry mean thickness `3.558764 mm` (SD `0.367833`, 2194 nonzero inner vertices); Lewiner geometry mean thickness `3.558218 mm` (SD `0.367726`, 2194 nonzero inner vertices); difference `-0.000546 mm`.

## Scope and limitations

- The harness compiles and invokes the selected `IsoSurfaceOnGrid` source directly. Minimal API stubs replace JIST progress/image containers and normal generation; the harness records the Java algorithm's actual returned vertices and triangle indices.
- Both methods consume identical float32 values. Per-volume SHA-256 digests are checked by the runner. The fixture's identity affine makes index coordinates equal to RAS millimetres.
- Vertex-nearest distances depend on tessellation density. Vertex-to-triangle distances here sample mesh vertices against the other mesh's triangles; they are not a continuous triangle-to-triangle Hausdorff bound.
- Although vertex coordinates and topology counts align closely, only the reported fraction of unordered faces match after mapping vertices; this indicates differing local connectivity/interpolation and is a relevant substitution risk for mesh-dependent measurements.
- Thickness uses the same Python closest-point operator on both extracted mesh sets. This quantifies geometry's effect on that measure; it is not a direct Java thickness implementation comparison.
- This single synthetic folded fixture does not establish clinical or end-to-end equivalence. Matching this case would not remove the need for broader shapes, adversarial topology cases, and source-version verification.

Full metric output and source/input provenance: [`comparison.json`](../../tests/data/synthetic/jist_surface_comparison/comparison.json). Fixture generator: [`make_phantom.py`](../../tests/data/synthetic/jist_surface_comparison/make_phantom.py).
