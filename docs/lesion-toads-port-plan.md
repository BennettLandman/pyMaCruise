# LesionTOADS and MaCRUISEplus port plan

**Status:** P39 is Started. P39a–P39zzh have validated bounded preparation, registration, membership-formula, spatial-neighborhood, fixed-input smoothing-composition, volume membership sweeps through OBJECT+LESION (including OPTIMIZED_LESION), final classification output, and thinning-frontier initialization, alongside BLACKHOLE prior/distance handling and combined OBJECT+LESION+OUTMASK, BLACKHOLE+OUTMASK, and LESION+BLACKHOLE+OUTMASK smoothing/blending/prior/distance/normalization,
atlas loading/sampling, registration derivatives, and selected registration
coefficient slices on synthetic/reference data, as well as LesionTOADS growing
queue exhaustion with pre-supplied memberships and P39t validates the
membership update's centroid-similarity setup. P39d/e validate bundled
compressed-NIfTI handling only; P39k reduces pre-sampled terms and does not
include atlas sampling, transform optimization, or segmentation. This
document tracks algorithm scope and source contracts; it does not claim a
full LesionTOADS port. P39zzn also matches repeated thinning with accepted and
rejected topology-driven label changes, and P39zzo matches the constrained
centroid updater with and without field correction; upstream membership
estimation remains unported. P39zzp also validates the selected weighted
centroid update for several blend factors, and P39zzq validates the selected
class/auxiliary variance update with its source branch behavior. P39zzr
validates selected diagonal-covariance channel weights and fallback behavior,
P39zzs validates pooled noise standard deviations, and P39zzt validates
generic centroid update modes. P39zzu validates source OUTLIER suppression in
the generic class-centroid accumulators, P39zzv validates PowerTable weights
for optimized CSF/GM/WM centroids, and P39zzw validates ordinary LESION and
BLACKHOLE centroid routing.
P39zzx validates its optimized LESION/BLACKHOLE variant and its separate
PowerTable factor source.
P39zzy combines the validated categories in one selected `computeCentroids()`
pass and validates cross-category linked updates.

## Selected workflow and source

The selected contract is the MaSeg-informed `maCRUISEplus_v1` path. MATLAB
`maCRUISEpluspipeline.m` invokes `run_LesionToads_pipeline.m`, which invokes
`run_LesionToads_jist.m` and the selected
`edu.jhu.ece.iacl.plugins.classification.MedicAlgorithmLesionToads` adapter.
The numerical driver is `AlgorithmLesionToads`, with
`LesionToadPreprocess` and `LesionToadSegmentation`; atlas handling and
registration are in `DemonToadDeformableAtlas` and
`DemonsToadsAtlasWarping` under the supplied `toads-cruise/src` tree.

The selected `run_LesionToads_jist.m` settings are: T1 MPRAGE plus FLAIR,
lesion modeling enabled, correction disabled, output type `dura removal
inputs`, maximum iterations 99, smoothing 0.2, energy tolerance 1e-4, atlas
prior 2.0, distance limits 3/25/2, and connectivity `(18,6)`. The adapter
sets rigid atlas alignment, four levels, 50 initial alignment iterations,
one main registration iteration, demons smoothing 1.0/scale 2.0, identity
normalization, lesion weighting, and prior centroid mode with smoothness 0.1.
The wrapper labels for the distance limits do not map one-to-one to the
constructor field names; verify downstream semantics before final parameter
claims. The sibling `CRUISEpluspipeline.m` requests correction and `cruise
inputs`; it is not the selected profile.

The legacy plus preprocessing also includes skull-mask handling, EM, and dura
removal. The selected Python contract requires a user-provided binary brain
mask on the T1 grid and excludes BET, SPECTRE, and FreeSurfer. The mask is
normalized to RAS with a warning when needed and is shared by T1/FLAIR
preparation after grid checks. EM and dura-removal behavior that feeds
LesionTOADS still requires source tracing and validation; this input decision
does not claim parity with those excluded skull-stripping stages.

## Completed substeps

### P39a: intensity and centroid initialization

`src/pymacruise/lesion_toads_preprocessing.py` ports negative-value clamping,
per-channel robust bounds (fraction 0.01, two refinement scales), division by
the robust span without subtracting the minimum, and lesion/black-hole
centroids `low/span + atlas_prior`. On the committed 9×8×7 T1/FLAIR phantom,
all 1,008 normalized values, both channels' bounds, and both centroid vectors
match the standalone Java harness exactly. Zero robust spans are rejected;
the Java source instead propagates NaN/infinity. This is an explicit
invalid-input deviation. See
`tests/data/synthetic/lesion_toads_preprocessing/README.md`.

### P39b: descriptor prior parsing

`src/pymacruise/lesion_atlas_priors.py` reads structure order, modality
intensity rows, lesion priors, and black-hole priors without loading image
arrays. A small image-free descriptor fixture checks ordering and malformed
inputs. The bundled `cruise-atlas-12obj-lesiontoads2012.txt` structure order
and selected T1 MPRAGE/FLAIR values match the source descriptor; lesion priors
are `[0.4, 1.2]` and black-hole priors `[0.1, 0.6]`. Missing selected profiles
are rejected instead of using the source's zero/label fallback.

### P39c: capped six-neighbor distance support

`src/pymacruise/lesion_distance.py::lesion_toads_distance_to_mask` ports only
the boolean-mask overload of `LesionToadSegmentation.fastDistance`. The
committed synthetic arrays compare every output voxel against a standalone
Java harness preserving the selected method's scan/cap behavior. This is a
six-connected voxel-step field, not a spacing-aware or Euclidean distance.
Objects are 0, reached background voxels have positive integer distance, and
unreached/empty-mask voxels remain -1. The source executes its initial scan
even with `max_distance=1`, so distance-2 voxels are assigned in that case;
the Python helper preserves this behavior. This validates only a support
primitive, not lesion membership or segmentation. See
`tests/data/synthetic/lesion_toads_distance/README.md`.

## Next source-faithful increments

1. Port transformed LesionTOADS prior sampling from the loaded `.nii.gz` shape
   maps. Match the selected class-specific bounds, rigid coordinate mapping,
   interpolation/default-value behavior, and preserve descriptor order. Keep
   this separate from registration and do not add an external atlas directory
   dependency.
2. Resolve and bundle the selected 18/6 critical-point lookup table as a
   regular Git file if the source algorithm needs it. It is absent from the
   selected `toads-cruise` source/atlas folders but present in the supplied
   VandyCRUISE source tree; record provenance and a checksum before use.
3. Port preprocessing, atlas alignment, lesion and black-hole membership
   initialization, relation/boundary maps, and the iterative
   topology-preserving competition/propagation stages in small observable
   increments. Reuse the validated intensity and descriptor helpers where
   their source contracts match.
4. Port the selected output products and the plus-specific preparation
   around the shared user brain mask, T1/FLAIR N4, EM, and dura stages.
5. Build a larger folded synthetic lesion phantom with a narrow sulcus,
   attached dura, and known lesion inclusion. Compare intermediate and final
   maps against the selected Java implementation; require per-label Dice
   ≥0.97 and record all surface/geometry limits when those outputs become
   available.
6. Assemble synthetic plus workflow outputs only after the base TOADS and
   required preparation stages can run independently. Keep actual MRI behind
   the synthetic readiness gate.

The selected algorithm is a separate LesionTOADS engine, so its preprocessing
and descriptor work can proceed while base P18 is incomplete. Full
MaCRUISEplus assembly remains dependent on shared input preparation and
verified tissue/dura products.

### Source-version and runtime caveats

Per owner direction, VandyCRUISE is the selected algorithm source for this
port. The exact deployed plugin JAR used by historical MaCRUISEplus runs has
not been identified, and companion `toads-cruise` copies differ in classes and
adapter details. Therefore the P39 slices validate only their named selected
VandyCRUISE source methods; no historical plus-workflow equivalence claim is
made until the deployed version is identified.

The selected VandyCRUISE `LesionToadSegmentation` does not call MIPAV
`AlgorithmMorphology3D` or `AlgorithmMorphology25D`. Those calls occur in the
alternate, older `toads-cruise/src/.../LesionToadSegmentation.java` copy
(SHA-256 `77016c000b8bbbfe55bfbe8d665ddd8be914c8ecbae8563ccce59aaa834bcfe0`),
while the selected VandyCRUISE file (SHA-256
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`) has no
MIPAV morphology calls. Because the deployed LesionTOADS version remains
unidentified, these alternate-version operations are not currently required
for the selected port. If source-version review makes them reachable, reuse
SciPy morphology primitives when their behavior matches; a Java
reimplementation is unnecessary when the library operation matches. Compare
the exact reachable call semantics, including dimensionality, structuring
element, border handling, dtype, and iteration behavior. Add a small adapter
or port only for behavior SciPy does not provide, and record any remaining
differences as explicit deviations.

### Completed substage: P39d topology NIfTI loading

`load_lesion_atlas_topology` resolves only the descriptor's topology `File`
entry and reads its bundled `.nii.gz` with the shared RAS-aware NIfTI loader.
It validates an in-atlas relative path and a 3D image, warns and reorients a
non-RAS input, and does not load any declared shape image. The bundled atlas
loads as an RAS-oriented `(172, 216, 164)` image. A 2×3×4 synthetic LPS image
checks voxel-axis reorientation and confirms an absent shape file is
irrelevant to this slice. P39zzf additionally checks an optional descriptor
dimension triplet against the RAS NIfTI shape and rejects malformed or
mismatched dimensions. Selected-source limitation: `LesionToadAtlas` reads
exactly `Nx*Ny*Nz` bytes from a headerless file, while the self-contained
repository bundles only `.nii.gz` topology data; this Python path therefore
uses NiBabel and does not claim to reproduce that raw-byte reader. Shape-array
loading, source index conversion, transformed sampling, and algorithmic
topology use remain open. Full-suite validation: 594 passed, 2 expected RAS
warnings, in 79.50 seconds on implementation commit
`f3ea4073b4933d0dc3a6a712f3bce59c9ec9763a`; see
[`docs/reports/full-pytest-2026-09-29-lesion-topology-geometry.md`](reports/full-pytest-2026-09-29-lesion-topology-geometry.md).

### Completed substage: P39e shape-prior NIfTI loading

`load_lesion_atlas_shape_priors` follows the descriptor's declared
`Structure`/`File` mapping and loads each of the 13 bundled shape maps as an
immutable float32 RAS array in declared structure order. It validates that
paths remain inside the atlas directory, that each image is 3D and matches the
descriptor dimensions, and that voxel values are finite. The shared image
loader warns/reorients non-RAS images; duplicate references such as
`OuterVentricular-WM` share the same loaded array as `Cerebrum-WM`. Actual
assets and synthetic malformed-path, missing-file, shape, dimension, value,
alias, and orientation cases are covered. This intentionally does not
implement the selected Java raw-float reader, atlas registration, or
transformed-prior interpolation; those remain separate steps.

### Completed substage: P39f transformed shape-prior sampling

`sample_lesion_shape_prior` follows the selected rigid branch of
`LesionToadAtlas.getRegularTransformedShape`. It reuses the source-compared
rigid coordinate mapping from `toads_shapes.py`, applies the LesionTOADS
per-class image bounds, then samples with the source `ImageFunctions`
constant-outside linear interpolation. Class index 0 uses default value 1;
other classes use 0. This path deliberately skips TOADS's closest-edge
interpolation and `shapeSlope` recovery. A headless selected-Java formula
harness covers identity/rigid mapping, anisotropic grid resolution, fractional
coordinates, per-class bounds, and low/high template out-of-range behavior.
The focused tests and full configured suite pass: 398 passed in 55.53 seconds
on tested implementation commit
`9a71937c62edb34c082f09c883624307d0236e9b`; see
[`docs/reports/full-pytest-2026-09-29-lesion-shape-sampling.md`](reports/full-pytest-2026-09-29-lesion-shape-sampling.md).
The harness does not instantiate the complete atlas class and this slice does
not perform registration, transformed-bound optimization, or lesion
segmentation. Provenance and fixture are in
`tests/data/synthetic/lesion_shape_sampling/README.md`.

### Completed substage: P39g transformed shape-prior derivatives

LesionTOADS registration differentiates each sampled prior in template
coordinates with `ImageFunctions.linearInterpolationXderivative`,
`Yderivative`, and `Zderivative`, then composes those values with the rigid
coordinate Jacobian. The selected helper already has a source-comparison
fixture in the TOADS registration tests; P39g makes that derivative available
through the LesionTOADS sampler while preserving its own prior-value and
input-bound behavior. Fractional and integer source-fixture derivatives,
including upper-boundary and out-of-range cases, match exactly. The full
suite passes 399 tests with two expected non-RAS reorientation warnings on
implementation commit `683ad5063549de8fcb56f40ec3cdabaab7505fab`; see
[`docs/reports/full-pytest-2026-09-29-lesion-derivatives.md`](reports/full-pytest-2026-09-29-lesion-derivatives.md).
This does not implement registration updates or segmentation.

### Completed substage: P39h one-voxel rigid-registration coefficients

`lesion_registration_voxel_contributions` composes the class-ordered
LesionTOADS sampler and derivative, the selected rigid image-to-template
Jacobian, and the shared source-compared coefficient accumulator. It follows
the selected `computeRegistrationCoefficients` class-0 and `priorT > 0`
gates. The 6³ synthetic formula case matches the Java excerpt exactly for
sampled priors, gradients, Jacobian, cost/norm, gradient, and diagonal
Hessian. Full suite: 401 passed in 56.22 seconds on implementation commit
`74e15e7a094a71dc14c1784d5f3bf0e1ec748d31`; see
[`docs/reports/full-pytest-2026-09-29-lesion-registration-voxel.md`](reports/full-pytest-2026-09-29-lesion-registration-voxel.md).
The fixture is documented in
`tests/data/synthetic/lesion_registration_voxel/README.md`. This composes a
single source-shaped voxel only: it does not implement scan-order volume
accumulation, energy normalization, transform optimization, atlas bounds
optimization, or lesion segmentation.

### Completed substage: P39i scaled registration-pyramid prior sampling

The LesionTOADS sampler now accepts registration pyramid scales and requires
explicit full-resolution input/template centers above scale 1. It reuses the
selected rigid coordinate map while sampling the caller-supplied coarse atlas
array with the lesion registration's constant-zero value and derivative
rules. On committed synthetic coarse maps, scales 2, 4, and 8 match the
selected-Java `ParametricTransform.imageToTemplate(..., scale)` and
`ImageFunctions` value/derivative formulas exactly, including anisotropic
resolutions, odd-dimension floor truncation, fractional samples, and both
constant-zero exterior directions. Full suite: 403 passed in 61.34 seconds
on implementation commit `bb35469f4849346bdf3e03fd466e7edcee7e5e9d`; see
[`docs/reports/full-pytest-2026-09-29-lesion-pyramid-scales.md`](reports/full-pytest-2026-09-29-lesion-pyramid-scales.md).
This is sampling support, not a complete atlas registration optimizer or
LesionTOADS segmentation port.

### Completed substage: P39j scaled single-voxel registration coefficients

`lesion_registration_voxel_contributions` now propagates pyramid scale through
prior interpolation and the selected rigid coordinate Jacobian before
composing the class-gated coefficient formula. The selected-Java harness
matches coordinates, priors, derivatives, Jacobians, cost/norm, gradient, and
diagonal Hessian exactly for scales 2/4/8, including class-0 exclusion,
multiple positive classes, zero prior, and out-of-template samples. The full
suite passed 404 tests in 52.36 seconds on implementation commit
`38db63cdc9ad85e0cb1b3097c4aec6c0af244842`; see
[`docs/reports/full-pytest-2026-09-29-lesion-scaled-registration-voxel.md`](reports/full-pytest-2026-09-29-lesion-scaled-registration-voxel.md).
This is a selected-formula single-voxel comparison only.

### Completed substage: P39k source-ordered coefficient accumulation

`accumulate_lesion_registration_coefficients` accepts pre-sampled
membership, prior, derivative, and Jacobian arrays, then reduces them in the
selected Java x/y/z/class order using float32 updates. Dense and offset-zero
step-two scans match the selected-method Java formula harness exactly for raw
cost/norm/gradient/Hessian and normalized measures. The reverse-order control
differs at float32 precision; zero norm preserves Java NaN division. The full
suite passed 406 tests in 60.56 seconds on implementation commit
`0740ef256011c30b225ca1d26ac72251b830d2bc`; see
[`docs/reports/full-pytest-2026-09-29-lesion-registration-volume.md`](reports/full-pytest-2026-09-29-lesion-registration-volume.md).
This isolates reduction semantics and does not sample the atlas or optimize
the transform.

### Completed substage: P39l fixed-transform scaled volume coefficients

`lesion_registration_volume_coefficients` samples each registered class map
and source spatial derivative at the selected coarse input voxels, computes
the scaled rigid coordinate Jacobian, and passes those terms through the P39k
source-order reducer. The scale-2 and scale-4 synthetic cases use odd
full-resolution image dimensions, anisotropic resolutions, a nonzero rigid
transform, four class maps, and step-two traversal. Raw and normalized
cost/norm/gradient/diagonal Hessian outputs match a harness using the selected
`ParametricTransform` and `ImageFunctions` methods exactly. The full suite
passed 407 tests in 52.44 seconds on implementation commit
`590bf3c4004fd7013ce6ae540bf8cddd4bf4be3f`; see
[`docs/reports/full-pytest-2026-09-29-lesion-scaled-registration-volume.md`](reports/full-pytest-2026-09-29-lesion-scaled-registration-volume.md).
This validates fixed-transform coefficient composition only, not transform
optimization or segmentation.

### Completed substage: P39m one scaled rigid-registration update

`register_lesion_shapes_scaled_step` uses the normalized P39l coefficient
vectors and delegates the diagonal update and strict-improvement lambda
search to the existing source-compared registration helper. It reevaluates
the scaled sampled-volume energy for each trial transform. Selected-Java
formula fixtures match all trial transforms/energies and final update fields
at pyramid scales 1 and 2; the latter backtracks eight times before accepting.
Full suite: 408 passed in 52.49 seconds on implementation commit
`f2216f08ee75a861eba86e6c4e7d5f6103a64c95`; see
[`docs/reports/full-pytest-2026-09-29-lesion-registration-step.md`](reports/full-pytest-2026-09-29-lesion-registration-step.md).
This validates one update, not repeated optimization or segmentation.

### Completed substage: P39n fixed-scale optimizer iteration loop

`optimize_lesion_shapes_scaled` carries the accepted transform and adaptive
lambda through repeated P39m steps while applying the selected `itMax`,
`minEdiff`, and `minLambda` loop conditions. The selected-Java synthetic
fixture matches every per-iteration line-search trace and final state. It
checks stop by lambda floor after six steps, an explicit three-step limit, and
a minimum-change threshold after one step. Full suite: 409 passed in 53.07
seconds on implementation commit
`76e6ec46a3940a3e5a584e725a9c751962e9e32e`; see
[`docs/reports/full-pytest-2026-09-29-lesion-registration-optimizer.md`](reports/full-pytest-2026-09-29-lesion-registration-optimizer.md).
This is a fixed-scale iteration slice, not pyramid orchestration or lesion
segmentation.

### Completed substage: P39o coarse-to-fine registration pyramid

`register_lesion_shapes_pyramid` composes the source-validated block-mean
subsampler and scale schedule with the P39n fixed-scale optimizer. A committed
odd 33³ multi-class phantom compares selected Java and Python outputs for
levels 1–4, including floor-truncated dimensions, per-level trial histories
and stop states, transform handoff, and the source lambda reset. The selected
method omits scale 1; levels 1 and 4 therefore have no optimization stage,
while level 4 runs scales 8→4→2. Zero-norm trial energies preserve Java NaN
behavior and serialize as JSON `null` in the fixture; those trials are
rejected by strict-increase logic. Full suite: 410 passed with two expected
RAS warnings in 54.14 seconds on implementation commit
`a194c16c3fe97aaacdb21797d8d85f0c79d19912`; see
[`docs/reports/full-pytest-2026-09-29-lesion-registration-pyramid.md`](reports/full-pytest-2026-09-29-lesion-registration-pyramid.md)
and [`tests/data/synthetic/lesion_registration_pyramid/README.md`](../tests/data/synthetic/lesion_registration_pyramid/README.md).
This is formula-level pyramid registration evidence only, not a full atlas
class comparison, lesion segmentation, or plus-workflow equivalence.

### Completed substage: P39p constructor-state initialization

`initialize_lesion_toads_segmentation_state` reproduces the selected
`LesionToadSegmentation` constructor's deterministic template-to-class maps,
source-order float32 class volumes and image centroids, and initial model
buffers. A committed 3×2×2 RAS NIfTI fixture compares a selected-constructor
formula harness for multi-class lesion and no-lesion cases. It checks absent
class `0/0` centroid NaNs, lesion-class selection, class similarity defaults,
signed-byte distance sentinels, membership/availability buffers, and identity
covariances. Six focused tests pass; the full configured suite reports 416
passed with two expected RAS warnings in 54.01 seconds on implementation
commit `c2b14be265412e7d011c02402a4f0f3bf4a78131`; see
[`docs/reports/full-pytest-2026-09-29-lesion-segmentation-state.md`](reports/full-pytest-2026-09-29-lesion-segmentation-state.md)
and [`tests/data/synthetic/lesion_segmentation_state/README.md`](../tests/data/synthetic/lesion_segmentation_state/README.md).
The harness is an excerpt of selected constructor formulas and does not
instantiate the full source class. No membership evolution, topology
propagation, final lesion segmentation, or plus workflow is validated.

### Completed substage: P39p constructor-state initialization

`initialize_lesion_toads_segmentation_state` reproduces the selected
`LesionToadSegmentation` constructor's deterministic template-to-class maps,
source-order float32 class volumes and image centroids, and initial model
buffers. A committed 3×2×2 RAS NIfTI fixture compares a selected-constructor
formula harness for multi-class lesion and no-lesion cases. It checks absent
class `0/0` centroid NaNs, lesion-class selection, class similarity defaults,
signed-byte distance sentinels, membership/availability buffers, and identity
covariances. Six focused tests pass; the full configured suite reports 416
passed with two expected RAS warnings in 54.01 seconds on implementation
commit `c2b14be265412e7d011c02402a4f0f3bf4a78131`; see
[`docs/reports/full-pytest-2026-09-29-lesion-segmentation-state.md`](reports/full-pytest-2026-09-29-lesion-segmentation-state.md)
and [`tests/data/synthetic/lesion_segmentation_state/README.md`](../tests/data/synthetic/lesion_segmentation_state/README.md).
The harness is an excerpt of selected constructor formulas and does not
instantiate the full source class. No membership evolution, topology
propagation, final lesion segmentation, or plus workflow is validated.

### Started substage: P39q growing-frontier initialization

The selected `initGrowingLabels` resets interior availability according to
skeleton emptiness, then seeds the min-heap with empty voxels that touch a
class on one of the six faces. Priorities compose `bestMembership` and
`maxMembershipRatio` with the voxel's existing ordering value. The lesion and
black-hole object types fold their separate membership maps into those scores
unless the optional max mode is enabled. The source's `bestMembership` code
compares the combined lesion score to its running maximum but stores the
ordinary membership as the next maximum; preserve this behavior when
reproducing the selected source. This slice stops at initialized queue and
availability state; queue events, topology checks, and the iterative
segmentation remain separate.

### Completed substage: P39q growing-frontier initialization

`initialize_lesion_toads_growing_frontier` reuses the source-compatible
min-heap and implements the selected `initGrowingLabels`, `bestMembership`,
and `maxMembershipRatio` behavior up to queue initialization. Four committed
7³ Java formula-reference cases match all availability bits and every popped
queue tuple for additive/max lesion scores, zero ties, and optimized-type
exclusions. The implementation preserves the source quirk that compares a
combined lesion score but stores ordinary membership in the running maximum.
Seven focused tests pass; the full suite reports 423 passed with two expected
RAS warnings in 53.99 seconds on implementation commit
`aee53ea15220aade46b83c546a2151f3fcb134e2`; see
[`docs/reports/full-pytest-2026-09-29-lesion-growing-frontier.md`](reports/full-pytest-2026-09-29-lesion-growing-frontier.md)
and [`tests/data/synthetic/lesion_growing_frontier/README.md`](../tests/data/synthetic/lesion_growing_frontier/README.md).
This does not validate queue transitions, topology evolution, or lesion
segmentation.

### Started substage: P39r one growing event

The source `propagateGrowing` pop reads the min-heap item, ignores it if the
skeleton center is already filled, otherwise applies the multi-label
18/6 homeomorphism check to the current hard segmentation. Accepted events
set the skeleton, segmentation, and ordering; update class volume using the
special lesion-class membership sum or maximum; clear center availability;
and add eligible six-face neighbors with lesion-aware priorities. Rejected
events restore availability at the event label. This next slice composes the
already validated local topology predicate and priority formulas and stops
after one queue pop; complete queue exhaustion/convergence is separate.

### Completed substage: P39r one growing event

`process_lesion_toads_growing_event` composes the selected 18/6 local
topology predicate and P39q lesion-aware scoring with one
`propagateGrowing` queue pop. Five committed 5³ reference cases match
accepted lesion updates in additive/max modes, non-lesion updates, topology
rejection, and stale-item handling. Comparisons include output label and
availability maps, exact float32 ordering/volume bits, inserted neighbors,
and the resulting min-heap pop sequence. Two focused tests pass; the full
suite reports 425 passed with two expected RAS warnings in 54.31 seconds on
implementation commit `b607fa74457e476da5337e4b7318a31dad2890f4`; see
[`docs/reports/full-pytest-2026-09-29-lesion-growing-event.md`](reports/full-pytest-2026-09-29-lesion-growing-event.md)
and [`tests/data/synthetic/lesion_growing_event/README.md`](../tests/data/synthetic/lesion_growing_event/README.md).
The reference composes a validated topology method body with selected
LesionTOADS event formulas; it is not a whole-class run or full segmentation
validation.

### Completed substage: P39s full growing-queue exhaustion

Compose the selected `initGrowingLabels` setup and `propagateGrowing` event
body until the source min-heap is empty. The synthetic phantom will contain a
folded class boundary and narrow bands so the run exercises source scan order,
lesion-aware priorities, accepted/rejected/stale events, and topology-gated
neighbor scheduling across many iterations. Compare the full event trace and
final segmentation, skeleton, availability, ordering, and class-volume state
against a selected-method Java harness. The committed RAS 9×9×9 folded label
phantom has three seeded class fronts and formula-generated float32
memberships, lesion/black-hole maps, orderings, and class volumes. Two Java
method-excerpt runs exhaust the source min-heap with additive and max lesion
scores (446 and 442 events). Python matches the full per-event traces,
including every pop priority and payload, topology result, and scheduled
neighbor, plus final segmentation, skeleton, availability, ordering and
class-volume bits, and empty residual queue. The harness uses the selected
18/6 LUT, source `Numerics.max`, and `BinaryHeap4D`; it is not a full
`LesionToadSegmentation` instantiation. The focused test passes; the full suite
reports 426 passed with two expected RAS warnings in 63.53 seconds on
implementation commit `46f1ef7ba47ee39f456c18d9a9bb6f9ecfac9dc8`. See
[`docs/reports/full-pytest-2026-09-29-lesion-growing-pass.md`](reports/full-pytest-2026-09-29-lesion-growing-pass.md)
and [`tests/data/synthetic/lesion_growing_pass/README.md`](../tests/data/synthetic/lesion_growing_pass/README.md).
Statistical membership estimation, competition/thinning, and final
LesionTOADS segmentation remain unvalidated.

### Completed substage: P39t centroid-similarity setup

The selected `computeMemberships` method first builds class-to-class and
lesion-to-class similarity coefficients from the current centroids, standard
deviations, channel weights, and `priorScale`. When T1 is available, a pair
involving `Sulcal-CSF` uses only the T1 channel; lesion similarity otherwise
uses all channels and excludes lesion, black-hole, and optimized-lesion
classes. The source's black-hole similarity calculation is commented out, so
this slice must preserve its existing state rather than invent a formula.
Four deterministic two-channel/four-class Java formula cases match the
class-similarity matrix and lesion-similarity vector bit-for-bit, including
the T1-specific Sulcal-CSF branch, the no-T1 fallback, lesion-type exclusions,
and zero-standard-deviation NaNs. The source-initialized black-hole similarity
vector is preserved unchanged because its only update in the selected method
is commented out. One focused test passes; the full suite reports 427 passed
with two expected RAS warnings in 63.35 seconds on implementation commit
`91fc287a4b1dffc1f4acc701e0c44beca2816d61`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-similarity.md`](reports/full-pytest-2026-09-29-lesion-toads-similarity.md)
and [`tests/data/synthetic/lesion_toads_similarity/README.md`](../tests/data/synthetic/lesion_toads_similarity/README.md).
This validates coefficient setup only; voxelwise membership estimation and
lesion segmentation remain open.

### Completed substage: P39u unfielded data energies

`compute_lesion_toads_unfielded_data_energies` reproduces the selected
`OBJECT`, `LESION`, and `BLACKHOLE` weighted squared-residual accumulators for
one voxel, before inversion, priors, or membership normalization. Nine
selected-Java cases compare all three output energy bits across three types
and three synthetic intensity vectors. The implementation preserves the
Java left-to-right float32 multiplication order. Two focused tests pass; the
full suite reports 429 passed with two expected RAS warnings in 54.25 seconds
on implementation commit `c69df55991a468279f1e31b213cb64783f910f24`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-data-terms.md`](reports/full-pytest-2026-09-29-lesion-toads-data-terms.md)
and [`tests/data/synthetic/lesion_toads_data_terms/README.md`](../tests/data/synthetic/lesion_toads_data_terms/README.md).
This does not cover field correction, Sulcal-CSF/MPRAGE or FLAIR-specific
PowerTable paths, smoothing, atlas priors, inversion, or normalized
memberships.

### Completed substage: P39v field-corrected data energies

Extend the checked unfielded `OBJECT`, `LESION`, and `BLACKHOLE` weighted
residual formulas to the selected `useField` branches by multiplying each
channel image by its corresponding bias field before subtracting its
centroid. Nine selected-Java cases match class, lesion, and black-hole energy
accumulator bits across three object types and three synthetic image/field
vectors. The focused tests report 2 passed; the full suite reports 431 passed
with two expected RAS warnings in 54.22 seconds on implementation commit
`c43ba3aaf7f38582dbd087bdb971f0d91bef2ebf`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-field-data-terms.md`](reports/full-pytest-2026-09-29-lesion-toads-field-data-terms.md)
and [`tests/data/synthetic/lesion_toads_field_data_terms/README.md`](../tests/data/synthetic/lesion_toads_field_data_terms/README.md).
This excludes special PowerTable handling, smoothing, atlas priors, inversion,
and membership normalization.

### Completed substage: P39w FLAIR lesion PowerTable energy

The selected lesion branch uses a special single-channel FLAIR residual when
`useLesionWeight` is enabled: square the float32 residual, pass it to the
cached lesion `PowerTable`, and add the returned double to the float32 lesion
energy accumulator (narrowing after addition). Six committed Java cases
cover cached-grid and out-of-range values with and without field correction;
class-energy and post-narrowing lesion-energy bits match exactly (maximum
absolute error 0.0). Two focused tests also check absent-FLAIR fallback and
the required exponent. The full suite reports 433 passed with two expected
RAS warnings in 54.03 seconds on implementation commit
`7ee767ca6c2c498115561ba5f57721e2ef3047e4`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-flair-weight.md`](reports/full-pytest-2026-09-29-lesion-toads-flair-weight.md)
and [`tests/data/synthetic/lesion_toads_flair_weight/README.md`](../tests/data/synthetic/lesion_toads_flair_weight/README.md).
This remains a data-term slice before smoothing, atlas priors, inversion, and
membership normalization.

### P39x Sulcal-CSF/MPRAGE data energy — source-correction follow-up

The selected source has distinct fielded and unfielded behavior for a generic
`OBJECT` named exactly `Sulcal-CSF`. In the unfielded branch, it accumulates
the special T1 squared residual first (or the CSF PowerTable lookup when
MPRAGE is present), then also executes the generic weighted all-channel object
energy loop: the apparent `else` before that loop is inside a comment. In the
fielded branch, the active `else` makes the special T1 term exclusive. If the
T1 channel is unavailable, the generic weighted multichannel term is used.

The corrected Python helper is pinned against a regenerated 10-case selected-
Java fixture, including a case that distinguishes source-order float32
accumulation from the reversed order. Focused regression: **10 passed in
0.75 seconds**, including data-term, field, FLAIR, and membership-composition
tests; `git diff --check` passed. The former nine-case description and full
suite result below refer to the pre-correction implementation. The old nine-
case evidence is superseded for source-branch correctness; a full-suite rerun
after the correction is pending. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-sulcal-csf.md`](reports/full-pytest-2026-09-29-lesion-toads-sulcal-csf.md)
and [`tests/data/synthetic/lesion_toads_sulcal_csf/README.md`](../tests/data/synthetic/lesion_toads_sulcal_csf/README.md).
This remains a data-term slice before smoothing, priors, inversion, and
membership normalization.

### Completed substage: P39y class membership conversion

The selected source converts an accumulated class energy with the strict
`num > ZERO` test: accepted energies become `1/num`, and all other values
become `INF`. A `DIRECT` relation stores this inverse unchanged; a `PRIOR`
relation multiplies by `shapeFactor(shape) = Numerics.min(1.0f, 2.0f*shape)`.
Eleven selected Java cases match output float32 bits for energies below, at,
and just above `ZERO`, direct and prior relations, shape factors below/at/above
saturation, and zero/negative energy sentinel behavior. The source
`Numerics.min` implementation determines tie behavior. Three focused tests
pass; the full suite reports 438 passed with two expected RAS warnings in
54.59 seconds on implementation commit
`6c92a3db7da65c96b2f01c0f06e791cd479b9813`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-membership-conversion.md`](reports/full-pytest-2026-09-29-lesion-toads-membership-conversion.md)
and [`tests/data/synthetic/lesion_toads_membership_conversion/README.md`](../tests/data/synthetic/lesion_toads_membership_conversion/README.md).
This remains one class-level conversion, not a full membership update.

### Completed substage: P39z lesion distance modulation

The selected source modulates lesion membership and `lesionRelation` using
valid ventricle distance first, or gray-matter distance only when the
ventricle branch is inactive and `GMflag` is true. A valid brainstem distance
then multiplies the current values independently. At zero ventricle/GM
distance, the numerator uses `spread`; brainstem uses its stored distance
directly, so zero clears the value. Preserve source maxima, `-1` sentinels,
branch precedence, and float32 operation order. This scalar slice does not
generate the distance maps or include black-hole membership.

Twelve selected-Java scalar cases match lesion membership and
`lesionRelation` float32 bits through ventricle/GM precedence and independent
brainstem modulation, including the special `spread` numerator at zero
ventricle/GM distance, `-1` sentinels, disabled maxima, and GM enable/disable.
The focused suite reports 2 passed; the full suite reports 440 passed with two
expected RAS warnings in 69.12 seconds on implementation commit
`c93ec3978753acf857a656c0d4792228d176c90d`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-distance-modulation.md`](reports/full-pytest-2026-09-29-lesion-toads-distance-modulation.md)
and [`tests/data/synthetic/lesion_toads_distance_modulation/README.md`](../tests/data/synthetic/lesion_toads_distance_modulation/README.md).

### Completed substage: P39za black-hole distance modulation

Port the selected `BLACKHOLE` branch's distance postprocessing. Apply the
selected ventricle-or-GM factor to `lesionMembership` and `lesionRelation`;
the source ventricle branch leaves `blackHoleMem` unchanged, while the GM
branch also scales it. A valid brainstem factor then independently scales
both auxiliary maps and `lesionRelation`. Preserve this source-specific
asymmetry, branch precedence, zero-distance `spread`, `-1` sentinels, and Java
float32 expression order. This reuses distances already computed elsewhere
and does not normalize voxel memberships.

Fourteen selected-Java cases match lesion membership, `blackHoleMem`, and
`lesionRelation` float32 bits. The cases confirm the source asymmetry: the
ventricle branch leaves `blackHoleMem` unchanged, while the GM and brainstem
branches scale it. Zero-distance `spread`, branch precedence, `-1` sentinels,
zero maxima, and GM control are covered. Focused suite: 4 passed; full suite:
442 passed with two expected RAS warnings in 81.37 seconds on implementation
commit `c5493e84aa72d03bd951194aad0a91eec0ee3e39`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-blackhole-distance-modulation.md`](reports/full-pytest-2026-09-29-lesion-toads-blackhole-distance-modulation.md)
and [`tests/data/synthetic/lesion_toads_blackhole_distance_modulation/README.md`](../tests/data/synthetic/lesion_toads_blackhole_distance_modulation/README.md).

### Completed substage: P39zb per-voxel membership normalization

Recreate the final one-voxel portion of `computeMemberships`: add lesion and
black-hole auxiliary membership values to the denominator in their source
branches, add each class membership in class order, normalize all class and
selected auxiliary values only when `den > 0`, and accumulate the local
energy from normalized class and auxiliary memberships. This slice begins
with assigned memberships and energies; it does not estimate their data
terms or spatial priors.

Eight selected-Java one-voxel cases match denominator, normalized class
memberships, lesion/black-hole/outlier auxiliary memberships, and local
energy float32 bits. They cover ordinary, LESION, BLACKHOLE, OUTLIER/OUTMASK,
nonpositive denominator, repeated auxiliary-map division, optimized-lesion,
and float32-sensitive sums. Focused suite: 5 passed; full suite: 447 passed
with two expected RAS warnings in 67.49 seconds on implementation commit
`4d0a25af19d146823402bea1febd5dfc60868695`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-voxel-membership-normalization.md`](reports/full-pytest-2026-09-29-lesion-toads-voxel-membership-normalization.md)
and [`tests/data/synthetic/lesion_toads_voxel_membership_normalization/README.md`](../tests/data/synthetic/lesion_toads_voxel_membership_normalization/README.md).

### Completed substage: P39zc selected lesion-membership integration

Compose the already checked per-class inverse/prior conversion, source
direct/prior auxiliary lesion relation initialization, lesion and black-hole
distance modulation, and final voxel denominator/normalization. Keep class
energies as synthetic precomputed inputs; smoothing and spatial-prior energy
assembly are not part of this integration slice.

Five selected-Java cases match class and auxiliary membership float32 bits,
`lesionRelation`, denominator, and local energy through class inversion,
direct/prior scaling, distance modulation, and normalization. They cover
mixed ordinary/lesion/black-hole classes, ventricle and GM paths, brainstem,
and zero denominator. Focused suite: 1 passed; full suite: 448 passed with
two expected RAS warnings in 61.93 seconds on implementation commit
`a55cc41d6c451f8d3698d1e80ee5128c9a909ee1`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-membership-integration.md`](reports/full-pytest-2026-09-29-lesion-toads-membership-integration.md)
and [`tests/data/synthetic/lesion_toads_membership_integration/README.md`](../tests/data/synthetic/lesion_toads_membership_integration/README.md).

### Completed substage: P39zd six-face smoothing-energy update

Port one source-ordered voxel update from the `smoothing > 0` branch of
`computeMemberships`. Sum, in the Java face/class order, all other classes'
neighbor membership squares and the lesion/black-hole auxiliary squares;
then apply the class-specific current-label self term and lesion/black-hole/
outmask energy cross terms. Preserve separate full and half-energy updates.
This does not iterate smoothing over a volume or update memberships.

Four selected-Java cases match full updates to the class, lesion,
black-hole, and outmask membership data terms, half-weighted class/lesion/
black-hole energy updates, and intermediate neighbor-square sums for
`OBJECT`, `LESION`, `BLACKHOLE`, and `OUTMASK`. The cases preserve face/class
loop order, current-class exclusion, lesion/black-hole cross-coupling, and
the distinct self-neighbor order. Focused suite: 1 passed; full suite: 449
passed with two expected RAS warnings in 52.76 seconds on implementation
commit `bba409dc9c892262067f8ceaeca06c5582d58e7b`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-smoothing-energy.md`](reports/full-pytest-2026-09-29-lesion-toads-smoothing-energy.md)
and [`tests/data/synthetic/lesion_toads_smoothing_energy/README.md`](../tests/data/synthetic/lesion_toads_smoothing_energy/README.md).

### Completed substage: P39ze smoothing-to-membership integration

Compose the P39zd full smoothing additions used for class/auxiliary
membership conversion and the separate half additions used for local energy
with the checked direct/prior conversions, distance branches, denominator,
and normalization. Keep the input a synthetic one-voxel stencil and keep
atlas-prior energy construction outside this slice.

Two selected-Java cases match normalized class and auxiliary membership
float32 bits, `lesionRelation`, denominator, and local energy after smoothing,
conversion, distance modulation, and normalization. The cases cover a
LESION/direct/ventricle/brainstem path and a BLACKHOLE/prior/GM-zero/brainstem
path. Focused suite: 1 passed; full suite: 450 passed with two expected RAS
warnings in 57.57 seconds on implementation commit
`eb0b7cc5d13a959608ad7acabc63675a494874cd`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-smoothing-membership-integration.md`](reports/full-pytest-2026-09-29-lesion-toads-smoothing-membership-integration.md)
and [`tests/data/synthetic/lesion_toads_smoothing_membership_integration/README.md`](../tests/data/synthetic/lesion_toads_smoothing_membership_integration/README.md).

### Completed substage: P39zf atlas-prior energy terms

Port the selected one-voxel `priorCoefficient` accumulation into class,
lesion, and black-hole membership data terms and their separate half-weighted
energy accumulators. Preserve source class-index exclusion and loop order,
lesion-class additions, and the special extra black-hole contribution when
the lesion class itself is a black-hole class. Inputs are already computed
shape/similarity vectors; building those fields and maps remains separate.

Six selected-Java one-voxel cases match class, lesion, and black-hole prior
sums; full data-term increments; and half-weighted energy increments for
ordinary, lesion, black-hole, optimized-lesion, and outmask targets. They
exercise absent/present lesion classes, a black-hole lesion class, the
`m != k` exclusion, zero priors/shapes, and float32-sensitive values.
Focused suite: 1 passed; full suite: 451 passed with two expected RAS
warnings in 57.63 seconds on implementation commit
`0a7fc6ca18fcc433568711839515f3e506201c10`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-prior-energy.md`](reports/full-pytest-2026-09-29-lesion-toads-prior-energy.md)
and [`tests/data/synthetic/lesion_toads_prior_energy/README.md`](../tests/data/synthetic/lesion_toads_prior_energy/README.md).

### Completed substage: P39zg smoothing/prior membership integration

Compose selected smoothing and prior additions to class and auxiliary data
terms and energy values before inverse/prior conversion, distance modulation,
denominator normalization, and local energy evaluation. Use committed
synthetic neighbor stencils and shape/similarity vectors; do not claim
similarity-field construction or volume-wide updates.

Two selected-Java cases match class and auxiliary membership float32 bits,
class energy, `lesionRelation`, denominator, and local energy through
smoothing, shape/similarity priors, conversion, distance modulation, and
normalization. They cover LESION/direct with ventricle/brainstem distances
and BLACKHOLE/prior with its extra black-hole lesion-class prior and
GM-zero/brainstem distances. Focused suite: 1 passed; full suite: 452 passed
with two expected RAS warnings in 57.47 seconds on implementation commit
`eb7e045a8927e3394a49e2f08c9352eae1f5cde0`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-smoothing-prior-membership.md`](reports/full-pytest-2026-09-29-lesion-toads-smoothing-prior-membership.md)
and [`tests/data/synthetic/lesion_toads_smoothing_prior_membership/README.md`](../tests/data/synthetic/lesion_toads_smoothing_prior_membership/README.md).

### Completed substage: P39zh OUTMASK mixture membership

Port the selected output assignment after the OUTMASK data and class terms
have been inverted. Direct relations blend class and outmask memberships by
the current outlier membership and update that outlier map; prior relations
apply the class shape factor to the mixed value and updated outlier map.
Keep later denominator normalization separate.

Twelve selected-Java scalar cases match class and outlier float32 output
bits across direct/prior relations, shape-factor range, outlier weights from
zero to one, and float32-sensitive inputs. Focused suite: 1 passed; full
suite: 453 passed with two expected RAS warnings in 57.68 seconds on
implementation commit `94a919564650fddd882e3cb30ba3dc6763e1a690`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-outmask-membership.md`](reports/full-pytest-2026-09-29-lesion-toads-outmask-membership.md)
and [`tests/data/synthetic/lesion_toads_outmask_membership/README.md`](../tests/data/synthetic/lesion_toads_outmask_membership/README.md).

### Completed substage: P39zi OUTMASK normalization integration

Compose the selected direct/prior class/outlier mixture with the shared
denominator and local-energy loop. Ensure the OUTMASK-written outlier map is
what normalization consumes; cover one and repeated OUTMASK class visits.
This uses already converted class/outmask terms and does not compute energies.

Three selected-Java cases match OUTMASK class values, outlier-map state,
denominator, normalized memberships, and local energy through direct/prior
mixture and shared normalization. They include zero and unit initial outlier
weights and two sequential OUTMASK classes. Focused suite: 1 passed; full
suite: 454 passed with two expected RAS warnings in 56.88 seconds on
implementation commit `899493c1c3999871b5c57f88a5df9a033f8701d9`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-outmask-normalization.md`](reports/full-pytest-2026-09-29-lesion-toads-outmask-normalization.md)
and [`tests/data/synthetic/lesion_toads_outmask_normalization/README.md`](../tests/data/synthetic/lesion_toads_outmask_normalization/README.md).

### Completed substage: P39zj lesion/black-hole plus OUTMASK voxel integration

Compose a lesion or black-hole auxiliary map and an OUTMASK/outlier map in
the same one-voxel denominator and normalization pass. Validate the
interactions with selected Java cases and preserve class ordering and
source-specific black-hole ventricle behavior. Memberships/energies begin as
preconverted scalar inputs.

Two selected-Java cases match class, lesion, black-hole, and outlier
membership float32 bits, relation, denominator, and local energy through
distance modulation, OUTMASK mixing, and shared normalization. They cover
LESION+OUTMASK with GM-zero/brainstem and BLACKHOLE+OUTMASK with ventricle
precedence (including the unchanged black-hole map) and a direct OUTMASK
relation. Focused suite: 1 passed; full suite: 455 passed with two expected
RAS warnings in 56.82 seconds on implementation commit
`edad578985013fe3a747a4fad42df11c22a73013`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-auxiliary-membership-integration.md`](reports/full-pytest-2026-09-29-lesion-toads-auxiliary-membership-integration.md)
and [`tests/data/synthetic/lesion_toads_auxiliary_membership_integration/README.md`](../tests/data/synthetic/lesion_toads_auxiliary_membership_integration/README.md).

### Completed substage: P39zk data-energy-to-membership integration

Compose selected unfielded class/lesion/black-hole residual energy
accumulation with the validated smoothing, shape/similarity prior,
membership conversion, distance, auxiliary-map, normalization, and local
energy kernels for a synthetic voxel. This is a one-voxel integration with
synthetic channel vectors and stencils; it does not run a complete volume.

Two selected-Java cases match computed unfielded class/lesion/black-hole
residual energies and final class and auxiliary memberships, denominator,
and local energy through smoothing, shape/similarity priors, conversion,
distance modulation, and normalization. They cover LESION and BLACKHOLE
targets with two synthetic channels. Focused suite: 1 passed; full suite:
456 passed with two expected RAS warnings in 64.21 seconds on implementation
commit `990407f49bf7ca271567f032a0d7d62a0bddfb06`. See
[`docs/reports/full-pytest-2026-09-29-lesion-toads-data-energy-membership.md`](reports/full-pytest-2026-09-29-lesion-toads-data-energy-membership.md)
and [`tests/data/synthetic/lesion_toads_data_energy_membership/README.md`](../tests/data/synthetic/lesion_toads_data_energy_membership/README.md).

### Started substage: P39zl field-corrected data-energy integration

Repeat the one-voxel data-energy-to-membership composition with source
`useField` enabled. Each synthetic channel sample is multiplied by its own
field value before residuals; all later smoothing, prior, distance, and
normalization steps stay in source order.

### Completed substage: P39zl field-corrected data-energy integration

Two selected-Java one-voxel cases compare field-corrected residual energies
and the resulting normalized memberships, denominator, and local energy for
LESION and BLACKHOLE branches. The synthetic fields have different nonunit
values per channel, and are applied before residual accumulation. All
compared float32 values match. The full suite reports 457 passed with two
expected RAS warnings in 56.82 seconds on implementation commit
`499de268c1ca4c1625d46d8b840e29ad673017c7`; see the
[full-suite report](reports/full-pytest-2026-09-29-lesion-toads-fielded-data-energy-membership.md)
and [fixture provenance](../tests/data/synthetic/lesion_toads_fielded_data_energy_membership/README.md).
This remains scalar integration using precomputed spatial/prior inputs.


### Completed substage: P39zm source-ordered six-face neighborhoods

`iter_lesion_toads_six_face_neighborhoods` streams class, lesion, and
black-hole face values for each interior voxel in the selected VandyCRUISE
`computeMemberships` x/y/z traversal. Each stencil is ordered X+, X-, Y+,
Y-, Z+, Z-. Voxels on the outer planes are not yielded, matching the Java
loop bounds. This helper never writes its inputs; the selected Java update
mutates membership arrays during its x/y/z traversal, so an eventual volume
update must explicitly preserve that in-place scan behavior. A Java formula harness and committed
5×4×6 coordinate-distinct float32 fixture compare every interior coordinate
and face value exactly. Focused tests: 2 passed. Full suite: 486 passed with
two expected RAS reorientation warnings in 80.05 seconds on implementation
commit `58a3d757e3a8eebb15d7f44cb72f44df416f27d`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-six-face-neighborhood.md`](reports/full-pytest-2026-09-29-lesion-toads-six-face-neighborhood.md)
and [`tests/data/synthetic/lesion_toads_six_face_neighborhood/README.md`](../tests/data/synthetic/lesion_toads_six_face_neighborhood/README.md).
This validates neighborhood extraction only, not the volume-wide membership
update or segmentation.


### Completed substage: P39zn streamed smoothing-energy composition

`iter_lesion_toads_volume_smoothing_updates` streams each six-face
neighborhood from P39zm through the source-compared scalar smoothing kernel
for every class, retaining the source interior x/y/z order and per-class
float32 operation order. The committed 5×4×6 coordinate-distinct Java
formula harness covers OBJECT, LESION, and BLACKHOLE terms with per-voxel
initial energy maps. Every output term at every interior coordinate matches
float32 bits; input volumes are unchanged. Focused suite: 2 passed. Full
suite: 488 passed with two expected RAS warnings in 86.61 seconds on
implementation commit `856e07b2246642b3dad298e43df52cdaa02eeb7b`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-volume-smoothing.md`](reports/full-pytest-2026-09-29-lesion-toads-volume-smoothing.md)
and [`tests/data/synthetic/lesion_toads_volume_smoothing/README.md`](../tests/data/synthetic/lesion_toads_volume_smoothing/README.md). This step
uses fixed memberships and does not implement in-place membership writes,
atlas priors, inversion, normalization, or segmentation equivalence.


### Completed substage: P39zo direct OBJECT in-place membership sweep

`run_lesion_toads_object_membership_sweep` applies source-ordered smoothing
terms from the current mutable six-face neighborhood, converts each ordinary
class with a direct relation, and normalizes all class values before moving
to the next interior voxel. A selected-Java formula harness matches every
interior post-normalization value bit-for-bit on a committed asymmetric
5×4×6, three-class fixture; the boundary planes remain unchanged and a
reverse-order scan produces a different result. Focused tests: 5 passed. Full
suite: 492 passed with two expected RAS warnings in 80.64 seconds on
implementation commit `5af89244adf65ffab25503e245c82184ee222d8c`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-object-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-object-sweep.md)
and [`tests/data/synthetic/lesion_toads_object_sweep/README.md`](../tests/data/synthetic/lesion_toads_object_sweep/README.md). Scope is limited
to positive shape support, direct relations, ordinary OBJECT classes, fixed
data terms, and smoothing; atlas priors and lesion/auxiliary types remain open.


### Completed substage: P39zp OBJECT shape-prior volume sweep

`run_lesion_toads_object_prior_membership_sweep` adds the selected pairwise
shape-prior energy to the ordinary OBJECT in-place sweep. It uses fixed class
priors and class-similarity coefficients plus per-voxel transformed shape
values, then follows the source direct/prior relation, shape gate, inversion,
and voxel normalization. The committed 5×4×6 Java formula fixture covers
mixed relations, coordinate-varying class shapes, and one zero-shape class;
all interior output membership bits match. Setting `priorCoefficient=0` and
all relations direct matches P39zo. Focused tests: 2 passed. Full suite: 494
passed with two expected RAS warnings in 79.30 seconds on implementation
commit `7b05f5c1da239cdfd0553e21e617da4b0d65441c`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-object-prior-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-object-prior-sweep.md)
and [`tests/data/synthetic/lesion_toads_object_prior_sweep/README.md`](../tests/data/synthetic/lesion_toads_object_prior_sweep/README.md). This
remains ordinary OBJECT-only evidence; lesion and auxiliary classes are open.


### Completed substage: P39zq OBJECT and LESION smoothing sweep

`run_lesion_toads_lesion_smoothing_sweep` extends the source-ordered
in-place update to ordinary OBJECT classes plus one LESION class. It updates
the separate lesion data term and auxiliary membership from the selected
six-face stencil, applies source shape gating and direct/prior conversion,
then normalizes each visited voxel. The selected-Java harness matches all
class and lesion membership bits on the committed 5×4×6 synthetic volume;
boundary planes are unchanged and a reverse-order control differs. Focused
tests: 3 passed. Full suite: 496 passed with two expected RAS warnings in
81.87 seconds on implementation commit
`db5f46743761cbe4feebbfb15343524b0787ae5a`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-lesion-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-lesion-sweep.md)
and [`tests/data/synthetic/lesion_toads_lesion_sweep/README.md`](../tests/data/synthetic/lesion_toads_lesion_sweep/README.md). Atlas priors,
distance modulation, black-hole and OUTMASK classes, and full segmentation
remain open.


### Completed substage: P39zr OBJECT+LESION shape-prior sweep

The P39zq volume sweep now optionally adds the selected pairwise class-prior
terms and lesion prior for one LESION class. It preserves the extra lesion
class contribution to ordinary class priors, lesionSimilarity accumulation,
source float32 order, shape gating, relation factors, in-place traversal,
and per-voxel normalization. On the committed 5×4×6 selected-Java formula
fixture, every interior class and lesion membership bit matches. The
zero-coefficient path matches P39zq. Focused tests: 4 passed. Full suite: 497
passed with two expected RAS warnings in 80.45 seconds on implementation
commit `e6eb60e8a962786343f498fb17e9b95f66673e43`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-lesion-prior-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-lesion-prior-sweep.md)
and [`tests/data/synthetic/lesion_toads_lesion_prior_sweep/README.md`](../tests/data/synthetic/lesion_toads_lesion_prior_sweep/README.md). Distance
modulation and other auxiliary class types remain open.


### Completed substage: P39zs lesion-distance volume routing

The OBJECT+LESION in-place sweep now optionally applies the selected
ventricle-first distance modulation, GM fallback, and independent brainstem
scaling after lesion membership conversion and before denominator
normalization. It retains distance-zero `spread`, signed-byte sentinels, and
lesion-relation scaling. The synthetic selected-Java harness matches every
interior class and lesion membership bit with GM flag enabled and disabled;
the no-distance path remains unchanged. Focused tests: 5 passed. Full suite:
498 passed with two expected RAS warnings in 91.20 seconds on implementation
commit `23eb94c768ad778e65004044f790c001124f43a7`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-lesion-distance-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-lesion-distance-sweep.md)
and [`tests/data/synthetic/lesion_toads_lesion_distance_sweep/README.md`](../tests/data/synthetic/lesion_toads_lesion_distance_sweep/README.md).
BLACKHOLE-specific distance interaction remains open.


### Completed substage: P39zt BLACKHOLE volume smoothing

`run_lesion_toads_blackhole_smoothing_sweep` adds one BLACKHOLE class to the
source in-place volume update, with separate lesion and black-hole data terms
and auxiliary memberships. It preserves the BLACKHOLE self-neighbor sums,
shape gate, relation conversion, combined auxiliary-before-class denominator
term, and per-voxel normalization. All interior outputs match a selected-Java
formula harness bit-for-bit on a committed asymmetric 5×4×6 fixture; outer
planes remain unchanged. Focused tests: 3 passed. Full suite: 500 passed with
two expected RAS warnings in 81.54 seconds on implementation commit
`089507f45ebc7ff659ee1dd3eb5c88fb6b0330fa`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-blackhole-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-blackhole-sweep.md)
and [`tests/data/synthetic/lesion_toads_blackhole_sweep/README.md`](../tests/data/synthetic/lesion_toads_blackhole_sweep/README.md). Black-hole
atlas priors and distance modulation remain open.


### Completed substage: P39zu BLACKHOLE shape-prior sweep

The BLACKHOLE volume sweep now optionally adds selected pairwise class,
lesion, and black-hole prior terms. It preserves the BLACKHOLE class's extra
lesion-class contribution to ordinary class priors, the separate lesion and
black-hole prior sums, float32 operation order, in-place traversal, relation
factors, shape gating, and denominator normalization. All class and auxiliary
membership bits match the selected-Java formula harness across the committed
5×4×6 volume; a zero coefficient matches P39zt. Focused tests: 3 passed. Full
suite: 501 passed with two expected RAS warnings in 80.53 seconds on
implementation commit `9fe0402f83748b3bd134fbfeee918f9ef4665f38`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-blackhole-prior-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-blackhole-prior-sweep.md)
and [`tests/data/synthetic/lesion_toads_blackhole_prior_sweep/README.md`](../tests/data/synthetic/lesion_toads_blackhole_prior_sweep/README.md). Distance
modulation remains open.


### Completed substage: P39zv BLACKHOLE distance volume routing

The BLACKHOLE sweep now routes ventricle, GM, and brainstem distance maps
through the selected auxiliary update. Ventricle precedence scales lesion
membership and lesion relation but leaves black-hole membership unchanged;
GM and brainstem scale both maps. The committed Java formula fixture compares
all interior class/auxiliary membership bits with GM enabled and disabled,
including zero-distance spread and independent brainstem scaling. Focused
tests: 4 passed. Full suite: 502 passed with two expected RAS warnings in
79.39 seconds on implementation commit
`c606d9c6c253b114d25581829d59dbec1e202a69`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-blackhole-distance-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-blackhole-distance-sweep.md)
and [`tests/data/synthetic/lesion_toads_blackhole_distance_sweep/README.md`](../tests/data/synthetic/lesion_toads_blackhole_distance_sweep/README.md).


### Completed substage: P39zw OUTMASK/outlier volume sweep

`run_lesion_toads_outmask_smoothing_sweep` updates ordinary OBJECT classes
plus one OUTMASK class in source x/y/z order. It preserves six-face
smoothing, per-class smoothing coefficients, direct/prior shape behavior,
outlier-weighted blending, outlier-map updates, shape gating, and shared
per-voxel normalization. The committed synthetic `5×4×6` Java formula
fixture covers zero, interior, and unit outlier weights; every interior class
and outlier-membership float32 bit matches. The harness preserves the Java
left-to-right `shapeFactor * outlier * numO` evaluation order. Focused tests:
5 passed. Full suite: 504 passed with two expected RAS warnings in 78.06
seconds on implementation commit `6895cf32d9af1ba989af2bd069997ca4bc62092f`;
see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-outmask-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-outmask-sweep.md)
and
[`tests/data/synthetic/lesion_toads_outmask_sweep/README.md`](../tests/data/synthetic/lesion_toads_outmask_sweep/README.md).
This is the selected OBJECT+OUTMASK formula branch, not a full
LesionTOADS-equivalence claim. Atlas similarity priors remain separate.


### Completed substage: P39zx OUTMASK atlas-prior volume sweep

The selected OBJECT+OUTMASK sweep now adds source-ordered class similarity
prior sums to the ordinary class energy and, for the OUTMASK class, to `numO`
before inversion. An asymmetric committed Java formula fixture matches every
interior class and outlier-membership float32 bit. A zero prior coefficient
reproduces P39zw exactly, and boundary voxels remain unchanged. Focused tests:
7 passed. Full suite: 506 passed with two expected RAS warnings in 77.67
seconds on implementation commit `ab80898125b2640f6f4bb118e4a6981448285ea5`;
see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-outmask-prior-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-outmask-prior-sweep.md)
and
[`tests/data/synthetic/lesion_toads_outmask_prior_sweep/README.md`](../tests/data/synthetic/lesion_toads_outmask_prior_sweep/README.md).
This validates the selected formula branch, not the full iterative
LesionTOADS engine.


### Started substage: P39zy LESION+OUTMASK volume composition

Compose the independently validated OBJECT+LESION and OBJECT+OUTMASK in-place
volume updates, including their selected smoothing and atlas-prior terms.
Preserve source ordering when the LESION auxiliary map is added to the
denominator before the class term, then the OUTMASK update blends and replaces
the outlier map before normalization. Distance modulation remains a separate
follow-on slice.


### Completed substage: P39zy combined LESION+OUTMASK volume sweep

The new `run_lesion_toads_lesion_outmask_smoothing_sweep` composes ordinary
OBJECT classes, one LESION class and its auxiliary map, and one OUTMASK class
with its outlier map in a single source-ordered pass. It combines six-face
smoothing, class and lesion atlas priors, direct/prior relations, OUTMASK
blending, shape gating, and shared denominator normalization. The committed
synthetic Java formula fixture matches every interior class, lesion, and
outlier float32 membership bit; boundaries are preserved, caller inputs remain
unchanged, and a zero prior coefficient reproduces the combined control.
Focused tests: 11 passed. Full suite: 509 passed with two expected RAS
warnings in 76.82 seconds on implementation commit
`c712c2ab983aa406d1983619b7c08ea4e446ba19`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-lesion-outmask-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-lesion-outmask-sweep.md)
and
[`tests/data/synthetic/lesion_toads_lesion_outmask_sweep/README.md`](../tests/data/synthetic/lesion_toads_lesion_outmask_sweep/README.md).
This validates one selected composed membership branch, not the full
iterative LesionTOADS engine. Distance modulation remains a separate step.


### Completed substage: P39zz combined lesion-distance routing

The combined LESION+OUTMASK sweep now routes the selected ventricle-first,
optional GM fallback, and independent brainstem distance rules through the
lesion auxiliary update, preserving the outlier blend and shared
normalization. Two committed Java formula fixtures cover direct and prior
lesion relations, ventricle precedence, GM fallback, zero-distance spread,
brainstem scaling, and varied outlier weights. Every interior class, lesion,
and outlier membership bit matches; zero distance maxima reproduce P39zy.
Focused tests: 10 passed. Full suite: 512 passed with two expected RAS
warnings in 77.36 seconds on implementation commit
`afed23eca5c3daec529339c7a00ec893a7dc8e21`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-lesion-outmask-distance-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-lesion-outmask-distance-sweep.md)
and
[`tests/data/synthetic/lesion_toads_lesion_outmask_distance_sweep/README.md`](../tests/data/synthetic/lesion_toads_lesion_outmask_distance_sweep/README.md).
This validates selected routing, not full LesionTOADS equivalence.


### Completed substage: P39zza BLACKHOLE+OUTMASK volume composition

The combined BLACKHOLE+OUTMASK volume adapter reproduces source-ordered
updates for either class ordering. Separate lesion/black-hole smoothing,
class and auxiliary atlas priors, the OUTMASK-only neighbor term, asymmetric
BLACKHOLE distance modulation, outlier blending, and shared normalization all
match two committed Java formula fixtures bit-for-bit for every interior
membership. Focused suite: 11 passed. Full suite: 515 passed with two
expected RAS warnings in 77.35 seconds on implementation commit
`1c89071df616d47135f2ee5c1690e36259815309`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-blackhole-outmask-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-blackhole-outmask-sweep.md)
and
[`tests/data/synthetic/lesion_toads_blackhole_outmask_sweep/README.md`](../tests/data/synthetic/lesion_toads_blackhole_outmask_sweep/README.md).
This remains selected formula validation, not whole-engine equivalence.


### Completed substage: P39zzb three-auxiliary volume composition

Compose LESION, BLACKHOLE, and OUTMASK classes in one source-ordered volume
sweep. Preserve distinct LESION and BLACKHOLE smoothing/prior equations,
distance-map behavior, OUTMASK blending, shared map updates, and denominator
normalization under multiple atlas class orders.

Two committed selected-Java synthetic references exercise both orders for the
last lesion-related class, which determines the source-selected `lesionClass`.
The fixtures cover all three auxiliary maps, smoothing and prior terms,
distance sentinels, zero/unit outlier weights, shape gating, shared map
normalization, and preserved boundaries. All interior class, lesion,
black-hole, and outlier float32 membership bits match. Focused suite: 13
passed. Full suite: 518 passed with two expected RAS orientation warnings in
87.97 seconds on implementation commit
`b44982bef1e9b87c34d597fa8c9fb9d008161908`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-three-auxiliary-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-three-auxiliary-sweep.md)
and
[`tests/data/synthetic/lesion_toads_three_auxiliary_sweep/README.md`](../tests/data/synthetic/lesion_toads_three_auxiliary_sweep/README.md).
This is selected formula validation, not whole-engine equivalence.


### Completed substage: P39zzc OPTIMIZED_LESION in three-auxiliary sweep

VandyCRUISE routes `OPTIMIZED_LESION` through the lesion membership, smoothing,
prior, inversion, and distance branches. Two selected-Java fixtures now test
this class before and after BLACKHOLE in the combined volume sweep, preserving
last-lesion-related-class selection and class iteration order. Both references
match all checked interior class and auxiliary float32 bits; duplicate LESION
and OPTIMIZED_LESION labels reject. Focused suite: 16 passed. Full suite: 521
passed with two expected RAS warnings in 84.10 seconds on implementation
commit `6c0835ede94ec411b0d22a82cbcc2c0045d6754d`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-optimized-lesion-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-optimized-lesion-sweep.md)
and the updated
[`tests/data/synthetic/lesion_toads_three_auxiliary_sweep/README.md`](../tests/data/synthetic/lesion_toads_three_auxiliary_sweep/README.md).
This validates the selected combined formula branch, not whole-engine
equivalence.


### Completed substage: P39zzd OPTIMIZED_LESION in LESION+OUTMASK sweep

The dedicated LESION+OUTMASK adapter now accepts `OPTIMIZED_LESION` and routes
its class, smoothing, prior, distance, conversion, and normalization behavior
through the shared lesion formulas. Two selected-Java synthetic fixtures cover
both class orders relative to OUTMASK and match every interior class, lesion,
and outlier membership bit. The fixtures also verify the source's
left-to-right prior outlier multiplication order; duplicate lesion aliases
reject. Focused suite: 15 passed. Full suite: 524 passed with two expected RAS
warnings in 84.97 seconds on implementation commit
`6fecb9cdf48573d6b61598c009c325a38249a699`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-optimized-lesion-outmask-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-optimized-lesion-outmask-sweep.md)
and the updated
[`tests/data/synthetic/lesion_toads_lesion_outmask_sweep/README.md`](../tests/data/synthetic/lesion_toads_lesion_outmask_sweep/README.md).
This validates the selected combined formula path, not whole-engine
equivalence.

### Completed substage: P39zzg final lesion/black-hole classification

`lesion_toads_final_classification` composes the selected loop's final
`updateClassificationFromLabels()` and `updateLesionClassifcation()` methods.
On the committed 3×2×2 RAS fixture, a selected-method Java excerpt and Python
match every final class label and every 0/1/2 lesion subtype, as well as the
changed-label and lesion-voxel counts. The cases exercise strict `>` tests,
lesion/black-hole ties, non-lesion class labels, and disabled lesion mode.
Focused tests: 4 passed. Full suite: 598 passed with 2 expected RAS warnings
in 87.88 seconds on implementation commit
`2b7515c415a6f4f65371c279d611ca4c50f4c5a3`; see
[`docs/reports/full-pytest-2026-09-29-lesion-final-classification.md`](reports/full-pytest-2026-09-29-lesion-final-classification.md)
and [`tests/data/synthetic/lesion_toads_final_classification/README.md`](../tests/data/synthetic/lesion_toads_final_classification/README.md).
This is final output mapping for supplied memberships only; upstream iteration
and plus-workflow equivalence remain open.


### Completed substage: P39zze OPTIMIZED_LESION in OBJECT+LESION sweep

The OBJECT+LESION volume adapter now accepts `OPTIMIZED_LESION` and uses the
source-selected lesion class formulas. Two Java references test middle and
last ordinary-class positions and match every checked interior class and
lesion float32 bit. A separate paired control shows exact output identity
between the aliases with atlas priors and distance modulation enabled.
Duplicate lesion classes reject. Focused suite: 8 passed. Full suite: 528
passed with two expected RAS warnings in 78.00 seconds on implementation
commit `312add6acbbefd7e570fa1a18c702432269bb71d`; see
[`docs/reports/full-pytest-2026-09-29-lesion-toads-optimized-lesion-object-sweep.md`](reports/full-pytest-2026-09-29-lesion-toads-optimized-lesion-object-sweep.md)
and the updated
[`tests/data/synthetic/lesion_toads_lesion_sweep/README.md`](../tests/data/synthetic/lesion_toads_lesion_sweep/README.md).
Java fixture coverage is the smoothing/normalization slice; the separate
operator references validate the prior and distance kernels.

### Completed substage: P39zzh LesionTOADS thinning frontier

`initialize_lesion_toads_thinning_frontier` ports the selected
`initThinningLabels()` setup: it resets class volumes, updates interior
skeleton and availability state, detects six-face segmentation boundaries,
and seeds the selected binary max heap using negative reciprocal
`maxMembershipRatio`. A source-method excerpt plus the copied selected Java
`BinaryHeap4D` match all resulting arrays, boundary counts, queue insertion
events, and queue removal order on the 4×4×4 synthetic RAS fixture for both
`useLesionMax` modes. This includes an interior non-boundary voxel, unchanged
outer planes, equal-priority heap ties, lesion/black-hole scores, and the
source `bestMembership()` quirk. Focused tests: 2 passed. Full suite: 600
passed with 2 expected RAS warnings in 80.34 seconds on implementation commit
`7820fc2d22ecee9ea531fc3c8c05aa44303fe1fe`; see
[`docs/reports/full-pytest-2026-09-29-lesion-thinning-init.md`](reports/full-pytest-2026-09-29-lesion-thinning-init.md)
and [`tests/data/synthetic/lesion_toads_thinning_initialization/README.md`](../tests/data/synthetic/lesion_toads_thinning_initialization/README.md).
This validates initialization only, not topology-gated changes during
`propagateThinning()` or full LesionTOADS segmentation.

### Completed substage: P39zzn topology-changing thinning pass

The selected thinning pass was checked with a second synthetic RAS NIfTI in
which memberships favor adjacent cyclic target labels across a folded,
mostly coherent tissue pattern. The Java reference invokes the selected
`CriticalPointLUT` predicate through the already compared 18/6 topology
adapter, then applies the selected heap, volume, and neighbor scheduling
semantics. Python matches all 706 heap events, candidate and topology
decisions (9 accepted, 69 rejected), 5 newly inserted neighbors, final
segmentation/skeleton/availability/ordering arrays, float32 class-volume
bits, and empty heap. Focused comparison: 1 passed. Full suite: 607 passed
with 2 expected RAS warnings in 89.92 seconds on implementation commit
`5c9b8d38833fe1834c96c5e8cc4fe4003076e570`; see
[`docs/reports/full-pytest-2026-09-30-lesion-topology-thinning.md`](reports/full-pytest-2026-09-30-lesion-topology-thinning.md)
and
[`tests/data/synthetic/lesion_toads_topology_thinning/README.md`](../tests/data/synthetic/lesion_toads_topology_thinning/README.md).
This is a fixed-input thinning-stage comparison; atlas sampling,
registration, upstream membership estimation, and full LesionTOADS remain
outside this validated slice.

### Completed substage: P39zzo constrained centroids

`compute_constrained_lesion_toads_centroids` ports the selected
`computeConstrainedCentroids()` method. It scans channels and classes in
source order and visits assigned voxels in x/y/z order, using membership
squared as the weight. With field correction enabled, each sample is
`field * image` before the weighted accumulation. It then applies the source
45% lower/upper bands using the neighboring centroids and the `-1e10` initial
sentinel. On the committed 5×4×3 RAS fixture, selected-method Java and Python
match raw float32 centroid bits with and without fields, including upper
clamps and an empty class whose zero mean is clamped to its lower band.
Focused tests: 2 passed. Full suite: 609 passed with 2 expected RAS warnings
in 79.97 seconds on implementation commit
`68dab53f1f851dcc46053fc4c0c317c6cff28f07`; see
[`docs/reports/full-pytest-2026-09-30-lesion-centroids.md`](reports/full-pytest-2026-09-30-lesion-centroids.md)
and
[`tests/data/synthetic/lesion_constrained_centroids/README.md`](../tests/data/synthetic/lesion_constrained_centroids/README.md).
This validates centroid updating only; other centroid categories,
variances/channel weights, and membership iteration remain open.

### Completed substage: P39zzp weighted centroids

`compute_weighted_lesion_toads_centroids` ports the selected
`computeWeightedCentroids(factor)` method. It accumulates
membership-squared-weighted image samples by current class in source order,
multiplies samples by their field values when enabled, and blends the mean
with the current centroid. The selected method sets an empty class to zero
without applying the blend, which differs from the constrained centroid
method and is preserved here. On the same committed 5×4×3 RAS phantom, Java
and Python match raw float32 bits for factors 0, 0.35, and 1 with and without
fields. Factor 1 preserves populated class values; the empty class remains
zero. Focused tests: 3 passed. Full suite: 610 passed with 2 expected RAS
warnings in 79.91 seconds on implementation commit
`8cb8ac17c57a54c4ece1e5400dcfb61500f41dce`; see
[`docs/reports/full-pytest-2026-09-30-lesion-weighted-centroids.md`](reports/full-pytest-2026-09-30-lesion-weighted-centroids.md)
and the shared
[`tests/data/synthetic/lesion_constrained_centroids/README.md`](../tests/data/synthetic/lesion_constrained_centroids/README.md).
This is a selected centroid method only; centroid-category assignment,
variance/channel-weight updates, membership updates, and complete segmentation
remain unvalidated.

### Completed substage: P39zzq class and auxiliary variances

`compute_lesion_toads_variances` ports the selected `computeVariances()`
method for class, lesion, and black-hole standard deviations. It accumulates
membership-squared residuals in source channel/class/x/y/z order and applies
field correction before residual calculation when enabled. On a 5×4×3
synthetic RAS phantom, the selected-method Java excerpt and Python match raw
float32 outputs for two channels, ordinary/LESION/BLACKHOLE classes, an empty
class, lesion and black-hole winner branches, ties, and zero auxiliary
memberships. The selected source black-hole winner branch assigns
`lmem = lesionMembership` before weighting the black-hole residual; the port
preserves and documents that behavior. Focused tests: 2 passed. Full suite:
612 passed with 2 expected RAS warnings in 78.04 seconds on implementation
commit `5f1948658030ec17407d1d53e3875d809b57cfda`; see
[`docs/reports/full-pytest-2026-09-30-lesion-variances.md`](reports/full-pytest-2026-09-30-lesion-variances.md)
and
[`tests/data/synthetic/lesion_variances/README.md`](../tests/data/synthetic/lesion_variances/README.md).
This validates the variance method for fixed maps only; membership iteration,
centroid-category assignment, and channel weights remain open.

### Completed substage: P39zzr diagonal channel weights

`compute_lesion_toads_channel_weights` ports selected
`computeChannelWeights()` behavior. For `Diagonal Covariance`, it accumulates
OBJECT, OUTLIER, LESION/OPTIMIZED_LESION, and BLACKHOLE terms, applies the
source `INF=1e10` when a channel numerator is nonpositive, then normalizes
weights over channels. It retains the unsquared linear lesion-related
residual sum and ignores MASK/OUTMASK types in this method. Other norm types
return uniform channel weights. The committed 5×4×3 RAS phantom matches raw
float32 Java values for diagonal mode with and without field correction,
uniform mode, and empty-contributor normalization. Focused tests: 2 passed.
Full suite: 614 passed with 2 expected RAS warnings in 90.66 seconds on
implementation commit `79ad6140610123725a1d758e4602aac8592a9b19`; see
[`docs/reports/full-pytest-2026-09-30-lesion-channel-weights.md`](reports/full-pytest-2026-09-30-lesion-channel-weights.md)
and
[`tests/data/synthetic/lesion_channel_weights/README.md`](../tests/data/synthetic/lesion_channel_weights/README.md).
This validates only channel weighting for fixed statistics and memberships.

### Completed substage: P39zzs pooled noise variance

`compute_lesion_toads_noise_standard_deviation` ports the selected
`computeNoiseVariance()` method. It scans hard classes 1..N−1, accumulates
membership-squared intensity residuals by channel in source order, and
broadcasts the pooled standard deviation to every class. Field multiplication
occurs before residual calculation when enabled. The existing 5×4×3 RAS
variance phantom matches selected-method Java bits with and without fields;
tests confirm background-class intensities do not affect the result and an
empty non-background pool yields all zeros. Focused tests: 2 passed. Full
suite: 616 passed with 2 expected RAS warnings in 113.03 seconds on
implementation commit `710085d2fe9407a27c77fcd641c17fa476c5ac33`; see
[`docs/reports/full-pytest-2026-09-30-lesion-noise-variance.md`](reports/full-pytest-2026-09-30-lesion-noise-variance.md)
and the shared
[`tests/data/synthetic/lesion_variances/README.md`](../tests/data/synthetic/lesion_variances/README.md).
This validates only the pooled estimator for fixed labels and centroids.

### Completed substage: P39zzt generic centroid modes

`compute_lesion_toads_basic_centroids` ports the generic class accumulator and
the selected `none`, `linked`, `damped`, `prior`, and fallback updates. It
retains source channel/class/x/y/z accumulation order, float32 narrowing,
optional field multiplication, and RAYLEIGH second-moment behavior for the
source modes that apply it. The 5×4×3 RAS variance phantom compares raw Java
float32 bits for two channels across all modes with and without fields. Its
classes cover ordinary OBJECT/MASK and RAYLEIGH categories plus an absent
class; link groups vary by channel, and the damped mode keeps the old centroid
for the empty class. LESION/BLACKHOLE and optimized LESION/BLACKHOLE types
reject until their specialized branches are ported. Focused tests: 2 passed. Full
suite: 618 passed with 2 expected RAS warnings in 123.29 seconds on
implementation commit `8275a1370ab76c8c5beab68518d65ab0372b6771`; see
[`docs/reports/full-pytest-2026-09-30-lesion-basic-centroids.md`](reports/full-pytest-2026-09-30-lesion-basic-centroids.md)
and the shared
[`tests/data/synthetic/lesion_variances/README.md`](../tests/data/synthetic/lesion_variances/README.md).
This validates generic centroid updates on supplied labels and memberships,
not specialized centroids or full LesionTOADS iteration.

### Completed substage: P39zzu OUTLIER centroid suppression

The generic centroid helper now applies the source OUTLIER adjustment
`mem = classMembership × (1 − outlierMembership)` before squaring and
accumulating. Its 5×4×3 RAS `.nii.gz` fixture contains two hard-labeled
OUTLIER classes and a map cycling through exact values 0, 0.25, 0.5, 0.75,
and 1. Raw float32 outputs match the selected Java method in all five centroid
modes, with and without field correction. Controls show OBJECT and OUTMASK
continue through the ordinary accumulator; absent class support, malformed
maps, and nonfinite maps are covered. Focused tests: 2 passed. Full suite:
620 passed with 2 expected RAS warnings in 110.09 seconds on implementation
commit `c87cbce7ad2023529cd64a4c018047815a34a05f`; see
[`docs/reports/full-pytest-2026-09-30-lesion-outlier-centroids.md`](reports/full-pytest-2026-09-30-lesion-outlier-centroids.md)
and the
[`tests/data/synthetic/lesion_outlier_centroids/README.md`](../tests/data/synthetic/lesion_outlier_centroids/README.md).
This validates OUTLIER centroid accumulation only, not outlier-map
estimation or full LesionTOADS iteration.

### Completed substage: P39zzv optimized tissue centroids

`compute_lesion_toads_basic_centroids` now applies the selected
`PowerTable.lookup(distance, exponent) / distance` weight for optimized GM
and WM classes and for optimized CSF only on the T1_MPRAGE channel. A CSF
class on FLAIR follows the generic accumulator. The method reuses the
source-compatible Python PowerTable implementation and preserves float32
residuals, factor narrowing, weighted accumulator order, field correction,
and each selected centroid update mode. The 5×4×3 RAS `.nii.gz` fixture
covers an exact zero residual, a squared residual below ZERO, cached table
lookup, direct out-of-range lookup, both CSF modality paths, an absent class,
and two channels. The harness compiles the preserved upstream Java
`PowerTable.java`; the source method and Python outputs match raw float32
centroid bits with and without fields for `none`, `linked`, `damped`, `prior`,
and fallback. Focused tests: 2 passed. Full suite: 622 passed with 2 expected
RAS warnings in 113.80 seconds on implementation commit
`92d6eb71e420986a404442c2b7a4dcdb6e3bce76`; see
[`docs/reports/full-pytest-2026-09-30-lesion-optimized-centroids.md`](reports/full-pytest-2026-09-30-lesion-optimized-centroids.md)
and
[`tests/data/synthetic/lesion_optimized_centroids/README.md`](../tests/data/synthetic/lesion_optimized_centroids/README.md).
This validates the optimized tissue centroid branch for supplied maps only.

### Completed substage: P39zzw LESION/BLACKHOLE centroid routing

`compute_lesion_toads_auxiliary_centroids` ports the ordinary LESION and
BLACKHOLE branches of selected `computeCentroids()`. It routes a hard-labeled
special-class sample by the maximum of class, lesion, and black-hole
membership; equality checks preserve class-first, lesion-second, black-hole
last tie precedence. Lesion-routed samples use the source lesion PowerTable
weight around the prior lesion centroid. The auxiliary lesion centroid then
blends 20% toward its weighted sample mean, while the black-hole centroid is
replaced by its direct weighted mean. The committed 5×4×3 RAS phantom covers
all three winners, three tie configurations, an exact zero lesion distance,
an absent class, and two channels with/without field correction. Selected Java
outputs match raw float32 bits for all five class-centroid modes, including
the source's hard-coded unlinked auxiliary updates. Focused tests: 2 passed.
Full suite: 624 passed with 2 expected RAS warnings in 155.83 seconds on
implementation commit `62d269854156c7fc45b2b35ab8f6c851e95b4e29`; see
[`docs/reports/full-pytest-2026-09-30-lesion-auxiliary-centroids.md`](reports/full-pytest-2026-09-30-lesion-auxiliary-centroids.md)
and
[`tests/data/synthetic/lesion_special_centroids/README.md`](../tests/data/synthetic/lesion_special_centroids/README.md).
This validates ordinary auxiliary-centroid routing for fixed memberships,
not optimized lesion categories or full iteration.

### Completed substage: P39zzx optimized LESION/BLACKHOLE centroids

`compute_lesion_toads_auxiliary_centroids` now supports
`OPTIMIZED_LESION` and `OPTIMIZED_BLACKHOLE` classes. It routes by the same
maximum-membership/tie rule, but computes a shared PowerTable factor from the
winning membership squared times field-corrected intensity, then uses it in
the selected class, lesion, or black-hole accumulator. This preserves the
selected source's unusual factor input rather than substituting a residual.
The RAS fixture adds a below-ZERO winning term and a PowerTable argument
outside the cached range. Selected Java outputs match raw float32 class and
auxiliary centroid bits for all five update modes with and without fields;
the routing coverage counts are 16 class, 7 lesion, and 7 black-hole samples,
with each of the three tie configurations present. Focused tests: 2 passed,
plus 2 ordinary auxiliary-centroid regression tests. Full suite: 626 passed
with 2 expected RAS warnings in 111.20 seconds on implementation commit
`9c785fa309bb38889736d6b494e0273e50182d87`; see
[`docs/reports/full-pytest-2026-09-30-lesion-optimized-auxiliary-centroids.md`](reports/full-pytest-2026-09-30-lesion-optimized-auxiliary-centroids.md)
and
[`tests/data/synthetic/lesion_optimized_auxiliary_centroids/README.md`](../tests/data/synthetic/lesion_optimized_auxiliary_centroids/README.md).
This validates optimized auxiliary centroid routing for fixed maps only.

### Completed substage: P39zzy mixed LesionTOADS centroid categories

`compute_lesion_toads_centroids` integrates the selected generic, OUTLIER,
RAYLEIGH, optimized CSF/GM/WM, ordinary LESION/BLACKHOLE, and optimized
LESION/BLACKHOLE accumulator branches in one shared class vector. It preserves
category-specific PowerTable inputs, field-corrected samples, class-first
auxiliary tie precedence, class/channel/voxel float32 accumulation, and
cross-category linked groups. The selected source applies the Rayleigh
square-root conversion only for centroid mode `none`; linked, damped, prior,
and fallback modes retain their direct expressions. The committed 6×5×4 RAS
fixture carries all categories, two channels and fields, both auxiliary maps,
and an intentionally absent class. Selected Java and Python outputs match
raw float32 bits for class, lesion, and black-hole centroids for all five
modes with and without fields. Focused centroid suites: 12 passed. Full
suite: 628 passed with 2 expected RAS warnings in 118.15 seconds on
implementation commit `fbecc45a357576d0c5a583d366d4dc5e1c3d6087`; see
[`docs/reports/full-pytest-2026-09-30-lesion-mixed-centroids.md`](reports/full-pytest-2026-09-30-lesion-mixed-centroids.md)
and [`tests/data/synthetic/lesion_mixed_centroids/README.md`](../tests/data/synthetic/lesion_mixed_centroids/README.md).
This validates a mixed fixed-map centroid update, not the enclosing
LesionTOADS membership iteration.

### Completed substages: P39zzzzk–P39zzzzm fixed-input composition

P39zzzzk applies the selected interior membership-map update to supplied
energy and spatial/prior inputs, including class, LESION, BLACKHOLE, and
OUTMASK maps. P39zzzzl separately reduces supplied per-class and auxiliary
energy maps into the selected objective. P39zzzzm builds raw unfielded
OBJECT/LESION/BLACKHOLE data-energy arrays over an asymmetric three-channel
RAS NIfTI, preserving x/y/z and class order plus shared auxiliary-slot
overwrites. Focused tests passed (1, 3, and 4 respectively). The full suite
reports 736 passed with 2 expected RAS warnings in 168.07 seconds on
implementation commit `a203b7ec912046b07017de34b63cb595a5f887af`; see
[`docs/reports/full-pytest-2026-09-30-toads-rayleigh-lesion-energy.md`](reports/full-pytest-2026-09-30-toads-rayleigh-lesion-energy.md).
These fixed-input adapters do not derive the full spatial/atlas energy state,
orchestrate all LesionTOADS passes, or establish workflow equivalence.
