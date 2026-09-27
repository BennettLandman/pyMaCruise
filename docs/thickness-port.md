# CRUISE thickness measurements

P31 follows two separate output branches from the V3 MATLAB workflow:

- `run_GetThickNessCRUISE.m` calls `GetThicknessCRUISE` (CVS r1.2). Its only
  enabled approximation is `Level Set Subtraction`, with iso-level 0,
  lambda 1, and step size 0.25. The selected `LevelSetThickness` solver (CVS
  r1.2) writes `abs(inner_levelset - outer_levelset)` at every voxel. Its
  historical mean, sample standard deviation, minimum, and maximum use only
  voxels where the inner and outer values are on opposite sides of zero.
- `run_GetThicknessByClosestPoints.m` calls `GetThicknessByClosestPoints`
  (CVS r1.3), which delegates to `ClosestPointThicknessForInner` (CVS r1.2,
  author Yuankai Huo). For each inner vertex, it measures the closest point
  on the outer triangles, then the return distance from that outer point to
  the inner triangles, and stores both values and their mean. When a central
  surface is supplied, it transfers the mean from the nearest vertex of the
  closest inner triangle and stores the central-to-inner triangle distance,
  nearest-inner-vertex distance, and vertex index.

`pymacruise.thickness.levelset_subtraction_thickness` and
`closest_point_thickness` implement these selected wrapper modes. Both return
physical millimetre values; image inputs are normalized to RAS with a warning
for non-RAS orientation, and image subtraction requires matching grids.
Closest points are computed against triangles rather than only mesh vertices.
A `scipy.spatial.cKDTree` over triangle centroids narrows candidates; the
global triangle-radius bound proves excluded triangles cannot be closer.
Thickness arrays are stored as VTK point-data fields and survive the local
ASCII VTK round-trip.

## Synthetic measurements

The committed fixture set `tests/data/synthetic/surfaces/thickness/` contains
0.5 mm RAS level sets for concentric spheres with analytic radii 5.5, 6.5, and
7.5 mm. Results measured with scikit-image triangulations:

| Measurement | Expected | Measured | Acceptance |
| --- | ---: | ---: | ---: |
| Level-set subtraction volume | 2.0 mm | 2.0 mm mean; 0.0 mm sample SD over 8,500 crossing voxels | Exact except float32 roundoff |
| Inner-surface symmetric closest-point mean | 2.0 mm | 1.9743 mm; 0.0387 mm sample SD | Absolute per-vertex error <0.15 mm |
| Central-to-inner closest-surface mean | 1.0 mm | 0.9961 mm | Absolute mean error <0.08 mm |
| Central transferred outer thickness mean | 2.0 mm | 1.9619 mm | Absolute mean error <0.12 mm |

The mesh sampling limits per-vertex distance accuracy; the volumetric
level-set subtraction remains exact for these parallel signed-distance
fields. Unit tests also verify non-RAS warning/reorientation, grid mismatch
rejection, unit outward normals, stable nearest-point correspondence indices,
and thickness-field VTK serialization. No actual imaging data are used.

## Scope and numerical interpretation

The wrappers' default level-set-subtraction thickness is not a Laplace or
streamline solution despite `lambda` and `stepSize` parameters being present;
those settings are accepted by the generic Java solver interface but do not
alter this selected method. Surface thickness is a symmetric closest-surface
distance, not normal-ray thickness. Rendering, visualization, and viewer-only
ASC export remain excluded.
