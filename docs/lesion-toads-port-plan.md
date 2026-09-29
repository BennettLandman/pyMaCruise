# LesionTOADS and MaCRUISEplus port plan

**Status:** P39 is Started. P39a–P39s have validated discrete preparation,
atlas loading/sampling, registration derivatives, and selected registration
coefficient slices on synthetic/reference data, as well as LesionTOADS growing
queue exhaustion with pre-supplied memberships and P39t validates the
membership update's centroid-similarity setup. P39d/e validate bundled
compressed-NIfTI handling only; P39k reduces pre-sampled terms and does not
include atlas sampling, transform optimization, or segmentation. This
document tracks algorithm scope and source contracts; it does not claim a
full LesionTOADS port.

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

### Started substage: P39d topology atlas loading

`load_lesion_atlas_topology` resolves only the descriptor's topology `File`
entry and reads its bundled `.nii.gz` with the shared RAS-aware NIfTI loader.
It validates an in-atlas relative path and a 3D image, warns and reorients a
non-RAS input, and does not load any declared shape image. The bundled atlas
loads as an RAS-oriented `(172, 216, 164)` image. A 2×3×4 synthetic LPS image
checks voxel-axis reorientation and confirms an absent shape file is
irrelevant to this slice. This intentionally does not emulate Java's
headerless `x + Nx*y + Nx*Ny*z` byte reader: it validates the repository's
compressed-NIfTI path only. Shape-array loading, source index conversion,
transformed sampling, and algorithmic topology use remain open.

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
