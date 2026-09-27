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

## Implemented TGDM kernel checkpoint

`src/pymacruise/tgdm.py` now provides a vectorized translation of the source
`GenericTGDM.updatePoint` local geometry and the `getEffectiveForce` methods
for the inner, central, and outer profiles. `tgdm_geometry` computes mirrored
boundary first differences, centered gradients, second/cross derivatives, and
the source mean/Gaussian curvature numerators. `tgdm_force` reproduces the
profile-specific normalization/transformation, pressure-weighted upwind
selection, outer-profile curvature cap/special case, and GVF advection term.
`tgdm_simple_point_618` reproduces the Java `SimplePointTest18`
geodesic-neighborhood encoding and the source six-connected interior plus
eighteen-connected exterior component checks used by the selected
`CONNECT_18_6` wrapper setting (including the constructor's mapping quirk).

The current focused tests check hand-calculated quadratic derivatives and
curvature numerators, flat-front pressure behavior, signed pressure upwinding,
signed advection upwinding, simple-point acceptance/rejection, and invalid
field rejection. The TGDM-specific suite has 8 tests. This validates isolated
source rules; it does **not** validate iterative TGDM results or numerical
equivalence to Java.

## Remaining source review and implementation

The Python kernel is only a building block. P29 still needs initialization and
signed-distance reinitialization, profile iteration schedules and time-step
switches, narrow-band construction/update ordering, landmine restarts,
bound/skeleton restrictions, tissue-force preparation for each surface, and
complete inner/central/outer array APIs.
After those are implemented, synthetic sphere/ellipsoid and nested-tissue
fixtures must check ordering, topology, stability, iteration accounting, and
all P08 metrics. Trusted Java reference outputs are still required for a
numerical equivalence claim.

The source has a connectivity inconsistency: `GenericTGDM` maps both
`CONNECT_6_18` and `CONNECT_18_6` to a six-connected foreground rule, while
`initNarrowBand` assigns the signs differently for those enum values. The
standalone wrapper requests `(18,6)`. Preserve this as an explicit source
question until the selected corrected V3 path/reference output settles the
intended behavior; do not silently substitute the independently used
topology-correction connectivity mapping.
