# Scientific level-set surface extraction

P30 follows MaCRUISE `run_LevelSetToSurf_jist.m` to VandyCRUISE
`GetMeshSurfaceFromLevelSet` CVS 1.4 and `IsoSurfaceOnGrid.solveOriginal`.
The selected wrapper passes connectivity 18, iso-level 0, and
`inclusive=false`. The Java source first subtracts the iso-level and moves
values in `(-0.1, 0]` to `-0.1`, and values in `(0, 0.1)` to `0.1`; this
near-zero rule is kept before extraction. The reference implementation uses
a connectivity-aware marching-cubes family. Python uses scikit-image's
Lewiner marching cubes implementation, which resolves ambiguous cubes and is
tested here for closed topology and physical geometry. Exact triangle and
vertex ordering is not part of this stage's contract.

`pymacruise.surfaces.extract_levelset_surface` accepts a 3D NIfTI image,
normalizes its voxel axes to RAS with a warning when needed, extracts the
zero set, then applies the complete voxel-to-world affine. The returned
`SurfaceMesh` uses RAS+ millimetres and triangular faces. VTK output is ASCII
legacy `POLYDATA` with an explicit RAS/mm title so scientific geometry can be
consumed without VTK/PyVista or rendering dependencies. Untagged legacy VTK
inputs need an explicit `input_frame`; LPS input is warned about, converted by
flipping R/A coordinate signs, and has face winding reversed to preserve
orientation.

The synthetic acceptance checks use analytic spheres and ellipsoids to
validate closed two-manifold edges, Euler characteristic 2, outward winding,
millimetre bounds, oblique affine handling, RAS normalization, LPS conversion,
and VTK round-trips. Their radius/ellipsoid errors quantify the triangulation
and voxel sampling error. On the committed 0.5 mm fixtures, the maximum
normalized ellipsoid-coordinate error is 0.0152 for the sphere and 0.0186 for
the ellipsoid; the oblique sphere has a maximum physical radial deviation of
0.0834 mm. This is scientific mesh data only; visualization
and rendering are out of scope (P47 discussion).

## Known scope

- Only ASCII legacy VTK triangular `POLYDATA` is read/written in this stage.
- The extractor preserves physical coordinates for oblique NIfTI affines and
  uses a marching-cubes triangulation. It does not duplicate JIST's exact
  table topology or mesh vertex ordering; compare the source mesh geometry on
  representative synthetic implicit fields before relying on subvoxel-level
  equivalence downstream.
- Non-RAS VTK files without frame metadata cannot be interpreted safely from
  coordinates alone. Callers must identify their frame explicitly.
