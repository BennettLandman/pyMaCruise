# CRUISE/TGDM port source map

P29 source tracing is underway. The MATLAB launch contract is
`run_MaCruise_jist.m` for the long CRUISE workflow and
`run_TGDMinvcrop_jist.m` for the standalone crop-space TGDM stage. The
selected VandyCRUISE implementation is `MedicAlgorithmNestedTGDM` plus
`GenericTGDM`; visualization and surface export are separate roadmap items.

## Standalone crop-space wrapper defaults

`run_TGDMinvcrop_jist.m` supplies membership input, tracking `Off`, all three
surfaces enabled, initial iso level 126.9, WM/GM iso levels 0.6/0.4, initial
smoothing curvature 1.0 for 6 iterations, inner curvature/pressure/external
forces 0.2/1/0 for 4 iterations, central forces 0.15/1.5/1 for 7 iterations,
outer forces 1/1/0, connectivity `(18,6)`, and generate iso-surfaces true. Its
script currently passes `IterCentralSurf` as the outer iteration count even
though the wrapper computes `IterOuterSurf`; this apparent wrapper defect must
be confirmed against the selected v3 workflow before porting.

The default standalone plugin itself differs: iso levels 0.5/0.6, initial
curvature 3.0, inner forces 1/1/0, outer curvature 0.25, and connectivity
`(18,6)`. P29 must preserve the wrapper's explicit values where supplied.

## Algorithm outline found so far

- If the initial image contains 0–255 memberships, the plugin computes
  `initIsoLevel - membership`, runs a fast-marching distance transform capped
  at 8, and smooths that initial surface before evolving it to the WM isosurface.
- It computes a signed tissue field from normalized WM and ACE-GM memberships.
  The central surface evolves from the inner result using that signed field and
  the optional GVF vector field.
- The outer surface evolves from the central result using gray matter, optional
  thinned ACE skeleton, and connectivity preservation.
- The plugin creates signed/level-set and surface products. Python porting for
  this roadmap stage covers level-set volumes; mesh extraction/output remains
  P30 and must report RAS+ millimetres.

## Remaining source review

Read and map `GenericTGDM` initial smoothing, inner/central/outer iteration
loops, narrow-band order, PDE time step, curvature/pressure/advection terms,
simple-point topology checks, speed-field normalization, brain-mask behavior,
and source fast-marching reinitialization. Add hand-calculated derivative and
force tests before full analytic sphere/ellipsoid evolution tests. No actual
imaging scan is part of this stage.
