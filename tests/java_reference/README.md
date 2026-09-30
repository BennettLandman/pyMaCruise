# Java algorithm reference harnesses

## ThinACE reference

`ThinACEReferenceHarness.java` runs the copied VandyCRUISE `ThinACE.java`
implementation on committed analytic masks. Regenerate the line, cube,
nested GM/WM, and symmetric sphere `.nii.gz` source/output fixtures with
`python scripts/generate_java_thinace_reference.py`. The harness uses minimal
test-only JIST image and voxel stubs under `stubs/`; it does not run a MIPAV or
JIST application. Source attribution and provenance are documented in the
fixture README and `license_info.md`.

The selected source applies `DistanceField.solve(4.0)`; the reference harness
uses the supplied MaCRUISE `DistanceFieldyh.java`, whose source header states
that it is identical to `DistanceField`. Temporary binary interchange files
and compiled classes are written under ignored `build/` and are not committed.

## Generalized GVF reference

`GVFReferenceHarness.java` runs the selected VandyCRUISE generalized GVF
implementation using source copies under `upstream/.../algorithms/gvf/`.
Regenerate the analytic ramp, asymmetric step, and Gaussian blob reference
vectors with `python scripts/generate_java_gvf_reference.py`. The selected
GVF, gradient, Poisson, multigrid transfer, and red/black solver sources retain
their upstream comments and CVS revision 1.1 metadata. Minimal JIST/framework
stubs exercise only the required in-memory image operations; they avoid
requiring MIPAV or a JIST application. Vector arrays retain the source
component order x/y/z, used by VandyCRUISE TGDM's corresponding x/y/z gradient
terms. Comparison metrics and the voxel-axis/physical-space convention are
documented in [`docs/gvf-port.md`](../../docs/gvf-port.md).

This test-only harness runs the selected Java MaCRUISE TGDM adapter on the
committed, deterministic 25³ synthetic phantom and writes small compressed
NIfTI (`.nii.gz`) golden outputs. It requires a JDK but does not run MIPAV or
JIST as applications. Regenerate the fixtures from the repository root with
`scripts/generate_java_tgdm_reference.sh`.

## Dura full-evolution reference

`DuraTGDMFullEvolutionReferenceHarness.java` compiles the selected VandyCRUISE
TGDM/profile source and runs it with the original JIST `DistanceField.class`.
Regenerate the small nested-surface and 25³ folded-cortex JSON checkpoints
with `scripts/generate_java_dura_original_jist_reference.sh`, setting
`JIST_TEST_JAR` and `JIST_LIB_DIR` to the supplied JIST build. This harness
uses the actual JIST image/distance classes with minimal headless framework
shims; it does not start MIPAV or the JIST application. The older
`scripts/generate_java_dura_tgdm_full_evolution_reference.sh` generator uses
the supplied `DistanceFieldyh` source as an auxiliary cross-check. On both
committed cases, those checkpoint arrays match the original JIST class
exactly. Detailed numerical metrics and the known voxel-level differences
between Python and Java are in `docs/toads-port-plan.md`.

## TOADS queue-step reference

`ToadsThinningQueueStepReferenceHarness.java` composes exactly one selected
`MultipleToadSegmentation.propagateThinning` event, using the actual copied
`CriticalPointLUT` and `BinaryHeap4D` classes. Seven synthetic 5³ cases compare
topology accepted/rejected and no-change paths, relation fallback, volume and
skeleton gates, output arrays, float32 class volumes, and the next heap order.
Regenerate with
`scripts/generate_java_toads_thinning_queue_step_reference.sh`. It is an event
slice only, not repeated thinning or complete TOADS validation. See the
fixture README and `docs/reports/toads-thinning-queue-step-2026-09-29.md`.

The `ToadsThinningQueueRunReferenceHarness.java` also replays the same selected
event under the original outer `val > -limit` loop, recording every pop and
state snapshot. Threshold fixtures disable neighbor insertion to isolate the
outer condition; a separate limit-100 case exercises dynamic rescheduling to
queue exhaustion. Its output is the synthetic fixture under
`tests/data/synthetic/toads_thinning_queue_run/`.


## Initialized TOADS thinning pass

`ToadsInitializedThinningPassReferenceHarness.java` composes source queue
initialization, minimum-volume calculation, the Java 18/6 LUT gate, dynamic
neighbor insertion, and outer-loop evolution on a 12³ synthetic phantom. It
matches all 1,000 seeded/popped items and final arrays to the Python pass.
Regenerate with
`scripts/generate_java_toads_initialized_thinning_pass_reference.sh`; the
fixture provenance and limits are in
`tests/data/synthetic/toads_initialized_thinning_pass/README.md`.

## Source used

`upstream/` contains source copies from the supplied checkouts, with their
original code comments and author/version metadata retained:

- MaCRUISE `MaGenericTGDM.java`, CVS revision 1.5, from the supplied
  `toads-cruise` tree.
- `DistanceFieldyh.java`, CVS revision 1.1, from the same tree. Its header
  describes it as identical to the omitted `DistanceField` class.
- VandyCRUISE `ProfileTGDM`, inner/central/outer/smoothing profiles, and the
  6-, 18-, and 26-connectivity implementations from the supplied source tree.
- JIST `BinaryMinHeap` and `Indexable`, including their original attribution
  comments, from the supplied JIST source tree.

The test-only `stubs/` provide the numerical subset of JIST image storage,
progress reporting, mesh return types, and unused framework methods needed to
compile those sources without bringing in the GUI/application dependency tree.
The reference uses the original TGDM solver, profile force equations, topology
rules, and distance-field algorithm. The harness maps the MaCRUISE adapter's
documented call sequence to those methods; it does not call mesh extraction,
which is tracked separately in P30.

These Java files are for regression comparison only and are not imported by
the Python package. Attribution and notice review is tracked in
[`license_info.md`](../../license_info.md); the supplied MIPAV/JIST notices are
preserved under `licenses/upstream/`.

## Phantom and outputs

The phantom is generated analytically in `TGDMReferenceHarness.java`: a radius
8 membership sphere, radius 5 white-matter sphere, synthetic gray matter,
zero GVF, and a planar thin skeleton. It has a 1 mm identity RAS affine. There
is no subject data. `tests/data/synthetic/tgdm/java_reference/` contains the
initial signed seed, smoothed seed, tissue-sign field, initial inner surface,
central surface, central shrink, final inner surface, and outer surface.


## TOADS growing-queue initialization

`ToadsGrowingQueueInitializationReferenceHarness.java` transcribes the selected
`initGrowingLabels` initialization against the test-only Java `BinaryHeap4D`.
Three synthetic cases verify class-specific six-face seed detection, float32
priority construction, availability changes, and full min-heap pop order.
Regenerate with
`scripts/generate_java_toads_growing_queue_initialization_reference.sh`; see
`tests/data/synthetic/toads_growing_queue/README.md`.

## TOADS topology-gated growing event

`ToadsGrowingQueueStepReferenceHarness.java` transcribes one queue event from
`MultipleToadSegmentation.propagateGrowing`, using the test-only Java
`BinaryHeap4D`, the selected 18/6 LUT, and the same homeomorphic predicate
harness used by the thinning references. Accepted, rejected, and stale-item
events record updated arrays, volume/order values, scheduled neighbors, and
the complete remaining heap. Regenerate with
`scripts/generate_java_toads_growing_queue_step_reference.sh`; synthetic
phantom details are in
`tests/data/synthetic/toads_growing_queue_step/README.md`.

## TOADS initialized growing pass

`ToadsGrowingPassReferenceHarness.java` composes source
`initGrowingLabels` with `propagateGrowing` until the Java min-heap is empty.
Its synthetic folded-band case records the full event trace and final arrays,
using the actual 18/6 LUT and `BinaryHeap4D`. Regenerate with
`scripts/generate_java_toads_growing_pass_reference.sh`; see
`tests/data/synthetic/toads_growing_pass/README.md` for the phantom and limits.
The same harness then refreshes relation map, boundary distances, and closest
labels on the post-growth segmentation, giving a Java reference for that
continuing-iteration call sequence.

## TOADS classification transfer from segmentation labels

`ToadsClassificationLabelsReferenceHarness.java` transcribes the selected
`updateClassificationFromLabels` loops over a small asymmetric categorical
volume and records the template-label lookup output plus changed-voxel count.
Regenerate with
`scripts/generate_java_toads_classification_labels_reference.sh`.

The Python-selected-method comparison and malformed-input checks are covered
by `tests/test_toads_classification_labels.py`; this reference does not imply
whole-algorithm equivalence.

## TOADS closest-label propagation

`ToadsClosestLabelsReferenceHarness.java` transcribes the selected
`computeClosestLabels` neighborhood seed and fixed-count diffusion loops.
It records all four signed-byte label slots on a folded asymmetric phantom at
three marching limits. Regenerate with
`scripts/generate_java_toads_closest_labels_reference.sh`.

## TOADS relation-list construction

`ToadsRelationListReferenceHarness.java` transcribes the selected
`buildRelationList` routine and records the compact relation map, relation
membership table, and per-class union rows for uniform, planar, and folded
synthetic label volumes. Regenerate with
`scripts/generate_java_toads_relation_list_reference.sh`.

## TOADS relation-map diffusion

`ToadsRelationDiffusionReferenceHarness.java` composes relation-list creation
with selected `computeRelations` exact-membership diffusion and records the
final map, input membership float bits, and loop-pass count. It reuses the
synthetic-label method excerpt in `ToadsRelationListReferenceHarness.java`.
Regenerate with
`scripts/generate_java_toads_relation_diffusion_reference.sh`.

## TOADS rigid registration update

`ToadsRigidRegistrationUpdateHarness.java` uses the supplied Java
`ParametricTransform`, `RotationMatrix`, and `ImageFunctions` with a compact
synthetic atlas/membership pair. It transcribes the selected exact-membership
coefficient loop and `registerGradientDescent` update for an accepted and an
11-trial backtracked case. Regenerate with
`scripts/generate_java_toads_rigid_registration_update_reference.sh`; see
`tests/data/synthetic/toads_rigid_registration_update/README.md` for scope and
limitations.

## TOADS boundary-distance update

`ToadsBoundaryDistancesReferenceHarness.java` transcribes
`computeBoundaryDistances`, including signed-byte interface seeding and the
source loop behavior that prevents a zero distance from expanding. Regenerate
with `scripts/generate_java_toads_boundary_distances_reference.sh`.

## TOADS relation-state refresh

`ToadsRelationStateRefreshReferenceHarness.java` calls the selected relation
diffusion, boundary-distance, and closest-label excerpts in wrapper order on
one folded synthetic input, carrying nonuniform float32 memberships and
prepopulated boundary bytes. It records every output and loop count.
Regenerate with
`scripts/generate_java_toads_relation_state_refresh_reference.sh`; fixture
scope is documented in
`tests/data/synthetic/toads_relation_state_refresh/README.md`.

## TOADS ordinary OBJECT intensity energy

`ToadsObjectIntensityEnergyReferenceHarness.java` records the selected
ordinary OBJECT squared-residual energy on a synthetic multimodal volume,
with optional field correction and PSF bias. Regenerate with
`scripts/generate_java_toads_object_intensity_energy_reference.sh`.

## TOADS ordinary OBJECT membership pass

`ToadsObjectMembershipPassReferenceHarness.java` composes the selected
ordinary-OBJECT intensity term, inverse shape prior, and ordered smoothing /
relation / reciprocal membership update on raw and field-corrected synthetic
inputs. Regenerate with
`scripts/generate_java_toads_object_membership_pass_reference.sh`.

## TOADS OBJECT/OUTLIER membership sweep

`ToadsObjectMaskOutlierSweepReferenceHarness.java` records a selected
source-order float32 membership sweep with ordinary OBJECT and OUTLIER
classes. It covers six-face smoothing, shape-relation factors, OUTLIER energy
mixing, normalization of both class membership and the scalar outlier map,
and objective accumulation on a deterministic 4 x 4 x 4 phantom. A paired
case also composes atlas inverse-prior energy and its half-weighted objective.
The generator checks the selected `MultipleToadSegmentation.java` source hash.
Regenerate
with `scripts/generate_java_toads_object_mask_outlier_sweep_reference.sh`;
the selected-method scope and limitations are listed in
`tests/data/synthetic/toads_object_mask_outlier_sweep/README.md`.

## TOADS OBJECT/MASK/OUTLIER membership pass

`ToadsObjectMaskOutlierMembershipPassReferenceHarness.java` combines the
selected OBJECT, MASK, and OUTLIER data branches with field correction, PSF
bias, six-face smoothing, atlas inverse prior, shape relation, normalization,
and objective accumulation on a 4 x 4 x 4 phantom. Regenerate with
`scripts/generate_java_toads_object_mask_outlier_membership_pass_reference.sh`;
the generator verifies the selected source hash. Fixture scope is documented
in `tests/data/synthetic/toads_object_mask_outlier_membership_pass/README.md`.

## TOADS MASK intensity energy

`ToadsMaskIntensityEnergyReferenceHarness.java` records the selected binary
`maskDistance(image, 0)` energy for positive, negative, exact-zero, and
signed-zero samples. Regenerate with
`scripts/generate_java_toads_mask_intensity_energy_reference.sh`.
