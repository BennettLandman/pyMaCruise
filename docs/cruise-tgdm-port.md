# CRUISE/TGDM port source map

P29 source tracing is underway. The selected base path is
`maCRUISEpipeline.m` → `run_maCruiseSurf_pipeline.m` →
`run_TGDMinvcrop_jist.m`, which launches the MaCRUISE adapter `MaTGDMcrop`.
That adapter uses the MaCRUISE `MaGenericTGDM` copy of the VandyCRUISE
`GenericTGDM` solver and profile classes. The standalone `MedicAlgorithmNestedTGDM`
contract is useful for source cross-checks but is not the selected wrapper.
Visualization and mesh extraction are separate roadmap items.

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

## Selected base MaCRUISE effective profile

The canonical `maCRUISEpipeline.m` sets `WMIsoLevel=0.5` and
`GVF_Lambda=1`; `run_TGDMinvcrop_jist.m` supplies `GMIsoLevel=0.4`, initial
iso 126.9, initial smoothing curvature/iterations 1.0/6, inner profile
0.2/1/0 for 4 wrapper iterations, central 0.15/1.5/1 for 7, and outer
1/1/0. The JIST command line accidentally connects `IterCentralSurf` to the
Java outer iteration argument instead of `IterOuterSurf`; both defaults are
7, so they coincide in the selected baseline. The wrapper request is
`CONNECT_18_6` and tracking is Off.

Unlike `MedicAlgorithmNestedTGDM`, the selected `MaTGDMcrop` then shrinks the
central level set for one inner-profile iteration with curvature input 0.15,
pressure -1, and no image force. It recomputes the delivered inner level set
from that shrunken central surface using the WM force with the source's
negated tissue field and pressure weight. The outer level set still starts
from the central surface. These adapter steps are included in
`tgdm_nested_arrays`.

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
field rejection. The TGDM-specific suite has 19 tests. It also checks zero-
force planar evolution, signed-pressure planar-front movement, bound handling,
wrapper force-field quantization, nested-stage routing, RAS warning/output
behavior, vector-component conversion, and a full nested run on the committed
25³ concentric-tissue phantom under reduced synthetic iteration settings. This
describes an earlier implementation snapshot; the authoritative complete-suite
pass count and tested commit are recorded in
[`reports/full-pytest-2026-09-27.md`](reports/full-pytest-2026-09-27.md). These
checks do **not** by themselves establish Java numerical equivalence.

`tgdm_evolve` now supplies the low-level source-style solver for an already
signed scalar level set. It applies the profile's narrow-band width, exact
wrapper-derived update schedule, bound clamps, local sign/topology checks,
landmine-triggered pass restarts, and signed-distance reinitialization. The
existing ACE fast-marching helper is reused; the supplied TOADS tree includes
`DistanceFieldyh.java`, explicitly described by its author as identical to
the `DistanceField` dependency used by TGDM. Source-output comparisons are
still required before asserting a match.

`tgdm_initial_membership_level_set`, the three pressure-field builders, and
`tgdm_nested_arrays` implement the selected wrapper adapters in array space:
membership-to-distance seeding and smoothing; the quantized 0..254 WM/ACE-GM
sign volume; the WM force; the central force; the GM/skeleton outer force; and
the initial inner → central → central shrink → final inner → outer sequence.
The wrapper passes all inputs on the same crop grid and applies its crop
metadata when writing/un-cropping; PyMaCRUISE will keep grid checks explicit.

`tgdm_nested_images` applies the shared NIfTI RAS normalization policy to all
scalar/vector inputs, emits warnings for axis reorientation, rejects differing
RAS grids instead of resampling, transforms GVF voxel-axis components with
axis flips/permutations, and returns RAS NIfTI images for every scalar output.
The source profile constructor accepts a `finalMaxSignDist` argument but does
not assign it to the field, so Java leaves every profile's final signed
distance cap at its field default of 20; the port preserves that behavior.

## Remaining source review and implementation

P29 still needs systematic synthetic nested-surface ordering/stability/
tolerance fixtures, larger-volume runtime and memory characterization, and Java
reference comparisons. Small default-schedule benchmarks on 25³–49³ phantoms
are recorded in [docs/tgdm-synthetic-benchmark.md](tgdm-synthetic-benchmark.md).
The selected adapter reads optional brain-mask and
ACE-skeleton inputs; these fields have been wired into bounds/topology
behavior, but need direct analytic tests.
Trusted Java reference outputs are required for a numerical equivalence claim.
The committed TGDM phantom and deterministic generator are under
`tests/data/synthetic/tgdm/`; `SHA256SUMS` covers all compressed NIfTI inputs.

The source has a connectivity inconsistency: `GenericTGDM` and its MaCRUISE
copy `MaGenericTGDM` map both
`CONNECT_6_18` and `CONNECT_18_6` to a six-connected foreground rule, while
`initNarrowBand` assigns the signs differently for those enum values. The
standalone wrapper requests `(18,6)`. Preserve this as an explicit source
question until the selected corrected V3 path/reference output settles the
intended behavior; do not silently substitute the independently used
topology-correction connectivity mapping.
