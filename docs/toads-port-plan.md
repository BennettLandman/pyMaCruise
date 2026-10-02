# TOADS, SPECTRE, and BET implementation plan

**Status:** P18 is in progress. P34 base orchestration is blocked until the
TOADS route below is implemented and validated. The user-supplied brain mask
decision removes BET/SPECTRE skull stripping from the runtime contract; it does
not establish algorithmic equivalence to either tool. The owner has selected
native T1 RAS space for the Python TOADS path. This is a deliberate working-
space deviation from historical V3, whose selected BET branch also ran
SPECTRE affine registration/resampling before TOADS.

## Selected implementation decisions

| Component | Decision | Runtime requirement | Validation evidence still required |
|---|---|---|---|
| TOADS tissue segmentation | Port the selected VandyCRUISE `MedicAlgorithmToads08` behavior and the MaCRUISE two-pass adapter contract into Python. Reuse NumPy/SciPy for array and image operations where they match source semantics; keep atlas alignment, iteration order, labels, masks, stopping rules, and wrapper settings source-derived. Do not substitute a generic classifier. | No MIPAV/JIST Java runtime in the finished package. Bundle all required atlas/configuration inputs under `data/atlases/` as regular Git files; image payloads stay `.nii.gz`. | Build a Java comparison harness around the supplied implementation. First compare intermediate/terminal fields on deterministic synthetic intensity, class-boundary, bias, and mask cases. Then compare final GM/WM-fill/WM-mask labels and geometry. Record voxel Dice per label, exact label agreement where deterministic, and all numeric deviations. |
| Dura removal between TOADS passes | Port the active `MedicAlgorithmRemoveDura09` stage and reproduce its selected placement between the two TOADS runs. Do not replace it with generic hole filling or morphology. | No separate external executable. | Use synthetic attached dura, disconnected outside tissue, thin bridges, and no-dura cases. Compare Java outputs and preserve exact masks/topology where the source is deterministic. |
| BET skull stripping | Exclude BET from the selected Python workflow. The user supplies the brain mask, which must be binary, nonempty, and on the T1 physical grid. | No FSL BET executable or Python BET wrapper. | Test binary/value validation, RAS reorientation warning, and exact shape/affine agreement against T1. This is an input-contract substitution, not a claim that the provided mask reproduces BET. |
| SPECTRE registration and skull-strip path | Exclude BET and SPECTRE from the Python runtime. The caller supplies a binary brain mask on the T1 grid, and the owner selected native T1 RAS as TOADS working space. Retain historical SPECTRE as a source/reference detail. Do not install or require FreeSurfer. | No FSL, BET, SPECTRE, or FreeSurfer runtime. | Historical V3 applies SPECTRE's 12-DOF affine registration with trilinear resampling before TOADS and later transfers labels back to the MaSeg grid. The selected Python path intentionally omits that transform/resampling. Record this in comparisons; use native-grid synthetic phantoms for stage validation and do not claim historical working-space equivalence. |

## Source inventory and available material

The active TOADS implementation is present in the supplied companion
`toads-cruise` tree, rather than the MIPAV checkout. The key code reviewed so
far is:

- `toads-cruise/src/edu/jhu/ece/iacl/plugins/classification/MedicAlgorithmToads08.java`:
  MIPAV/JIST wrapper, selected output type, active direct-copy image input,
  and wrapper parameter defaults/forwarding. Its isotropic WSINC resampling
  call is commented out and is not part of the selected execution path.
- `toads-cruise/src/edu/jhmi/rad/medic/algorithms/AlgorithmMultipleToads.java`:
  1,282-line workflow/algorithm driver, including input preprocessing,
  inhomogeneity correction, atlas alignment, segmentation and result products.
- `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadDeformableAtlas.java`:
  3,184-line deformable atlas implementation.
- `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadPreprocess.java`:
  960-line preprocessing implementation.
- `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java`:
  3,324-line classification/segmentation implementation.
- `VandyCRUISE/src/edu/jhu/ece/iacl/plugins/segmentation/skull_strip/MedicAlgorithmRemoveDura09.java`:
  selected inter-pass dura-removal wrapper and parameters.
- `data/atlases/cruise-atlas-10obj/` and `data/atlases/lesionToads-atlas-2012/`:
  bundled atlas descriptions/images already staged in PyMaCRUISE.

The listed implementations use MIPAV image classes as data containers and
infrastructure. The selected Python path must replace those containers with
RAS NIfTI/NumPy/SciPy types and map each numeric algorithm operation; the MIPAV
checkout is not a runtime prerequisite. No independent Python implementation
of this selected source/version has yet been identified. The port must not
shrink the work to a generic intensity classifier.

## Bundled TOADS atlas geometry discrepancy

The selected root descriptor declares a 172×216×164 topology grid at
0.9375 mm. Its source MIPAV XML and the bundled topology `.nii.gz` instead
declare 1.0 mm. This is now resolved as a source contract: the selected
`MultipleToadDeformableAtlas.loadAtlas` parser reads the topology and shape
resolutions from the text descriptor (source SHA-256
`df2f2ace677f2222e94c946730adf485d4401a76f7a9d2b302c9e5c8b652596c`, topology
resolution assignments at lines 986–987 and per-class shape resolution at
1061), then `loadTemplateImage` and `loadShapeImage` read headerless byte/float
payloads using descriptor dimensions (lines 818–859, 990, 1077). The Java
algorithm therefore uses descriptor spacing, not the NIfTI affine or source
MIPAV XML. `load_toads_atlas` retains both metadata values and emits a visible
warning for the topology discrepancy; it never silently resamples the atlas
or rewrites its NIfTI header. This explains the selected loader behavior but
does not establish full atlas geometry alignment or segmentation parity.
The root descriptor's filenames were changed to `.nii.gz` for the Python
repository's storage contract; the Java source reader expects raw payloads.
The Python source-index array loader is implemented and exact against the 12
available raw source maps. Thus no source raw files need to be bundled: runtime
loading uses only repository NIfTIs and the committed conversion audit. The
descriptor-versus-header spacing contract is documented in
[the P18ac source audit](reports/toads-atlas-spacing-source-audit-2026-09-28.md).
Physical atlas-to-input alignment remains unvalidated and is separate from the
confirmed source loader rule.

### Completed substage: fixed-transform shape-prior sampling

`src/pymacruise/toads_shapes.py` ports the selected rigid transform mapping,
closest-edge trilinear interpolation, class-specific input-bounds defaults,
and `shapeSlope` recovery/clipping around
`MultipleToadDeformableAtlas.getRegularTransformedShape`. The six-value rigid
parameter vector uses the source inverse-Rodrigues transform and millimetre
translations; input arrays use RAS-oriented x/y/z axes, and the template
spacing passed to the map comes from the atlas descriptor (0.9375 mm in the
selected atlas). Eight synthetic cases match the copied Java
`ParametricTransform`, `ImageFunctions`, and `Numerics` APIs exactly at
float32 precision, including mapped template coordinates. This fixture
transcribes the private atlas method's class-bounds and slope wrapper; it does
not instantiate the atlas class, validate class-bound optimization, perform
registration, or resolve the converted NIfTI world origin. Registration
optimization and physical-world geometry reconciliation remain separate work.

The selected TOADS segmentation also requires `critical186LUT.raw.gz` for
18/6 connectivity. The 936 KB compressed asset is bundled at
`src/pymacruise/data/critical186LUT.raw.gz`; its compressed SHA-256 is
`cbf1e339d78a0a2c673c939ee2c8a30ee3a35bfd9603d07ed2d48864bb0e88be`. Its
67,108,864-byte payload is indexed by a
26-bit x-major/y/z neighborhood key with the center omitted. This table is
not interchangeable with P25's current 6/18 predicate: the two already differ
for the single-corner key. `critical186_is_regular` now applies the bundled
table to a single binary neighborhood. Sixty-eight committed synthetic cases
match the copied VandyCRUISE `CriticalPointLUT` key and lookup methods exactly
(see [the fixture README](../tests/data/synthetic/toads_critical_lut/README.md)
and P18f/P18h in `plan.md`). This local LUT query does not port
`MultipleToadSegmentation.isHomeomorphicSegmentation`, which also checks old
and new objects and neighboring label relations, pairs, and triples.

### Completed substage: 18/6 multi-label homeomorphic-change gate

`src/pymacruise/toads_topology.py` ports the selected
`MultipleToadSegmentation.isHomeomorphicSegmentation` local accept/reject
logic with `checkComposed == false`: unchanged center labels accept
immediately; otherwise the 18/6 LUT checks the proposed object, the old
object, each distinct neighboring-label relation, and every unordered pair
of neighboring-label relations for both old and new labels. The input is a
complete x/y/z-ordered 3×3×3 short-label neighborhood. Callers must provide
an interior voxel or padding, matching the source's unchecked 26-neighbor
access.

All 130 committed synthetic cases match the copied selected Java method-body
excerpt against the actual selected `CriticalPointLUT` source and bundled
18/6 table. This is a local topology gate only: it does not port TOADS
membership evolution, priority/scheduling, multi-voxel update ordering,
`checkComposed` variants, border policy, or full segmentation behavior. See
the fixture notes at
`tests/data/synthetic/toads_homeomorphic_change/README.md` and the harness in
`tests/java_reference/ToadsHomeomorphicChangeReferenceHarness.java`.

### Completed substage: transformed shape support bounds

`compute_toad_shape_bounds` in `src/pymacruise/toads_shapes.py` ports the
selected rigid `precomputeImageToTemplateMatrix` path and the atlas wrapper's
voxel scan for positive interpolated shape support. Its five synthetic
6×5×4 cases match the copied Java transform/interpolation APIs and the
`computeTransformedShapeBoundingBox` source-body excerpt at every positive
support coordinate and bound. They cover anisotropic voxel sizes, rigid
rotation/translation, full support touching every image face, and an empty
support atlas prior.

The source stores the greatest support index as its `max` value without
adding one, then tests that value with an exclusive `>= max` condition during
shape sampling. PyMaCRUISE preserves and tests this historical behavior,
including the source empty-support sentinel (`min=input dimensions`,
`max=0`). Changing it would require a separate source-faithfulness decision.
This is only the fixed-transform support scan; registration optimization,
iterative bound updates, and physical atlas-world origins remain open. See
`tests/data/synthetic/toads_shape_bounds/README.md` for the fixture boundary.
The active exact-membership rigid objective, its update rule, and the
ANTsPy/SciPy substitution assessment are tracked in
[toads-atlas-registration.md](toads-atlas-registration.md).

## Historical V3 wrapper contract and selected Python deviation

The historical `run_maCRUISE_v3.m` default selects BET. Its `run_Bet.m` branch
runs BET, then SPECTRE affine registration of that stripped image to the
OASIS atlas using correlation ratio and trilinear interpolation. TOADS pass 1
consumes the SPECTRE-transformed image; dura removal's `inOriginal` is also
that stripped/registered image; pass 2 consumes the dura-stripped result; and
the three TOADS products are transferred to the MaSeg reference grid. The
SPECTRE step therefore changes TOADS working space and image sampling, in
addition to skull-mask preparation.

The selected Python contract is explicit: caller-provided binary mask,
native T1 RAS working space, and no BET/SPECTRE resampling. N4 T1 values are
zeroed outside that mask before TOADS, and the same masked native-space image
is supplied as `inOriginal` to dura removal. This preserves the selected
algorithm/atlas/parameter/output contracts while intentionally changing the
historical working space and trilinear interpolation effects. Routine
TOADS-output transfer-back is identity on this path; the generic nonidentity
transfer utility remains tested separately. P37/P43 must report the working-
space difference rather than attributing it to numerical noise.

`run_Toads_jist.m` passes `-inAtlas2 2.0`, `-inSmooting 0.2`,
`-inMaximum 1e-4`, `-inMaximum2 50`, rigid atlas alignment, and `(18,6)`
connectivity. In the selected `toads-cruise` constructor these map to
`atlasCoefficient=2.0`, `atlasScale=0.2`, `smoothing=0.2`,
`maxDistance=1e-4`, and `iterations=50`; `2.0` is not the atlas scale.
Additional constructor values are `marchingLimit=3.0`, `volumeLimit=0.5`,
`levels=4`, `firstIter=50`, `topIter=5`, `lastIter=0`,
`correctInhomogeneity=true`, `Chebyshev`, polynomial degree `3`, and kernel
size `30.0`. The wrapper also fixes `inOutput2=false` and
`inOutput3=false`. These values were traced through `MedicAlgorithmToads08`
into the selected companion `AlgorithmMultipleToads` constructor. `run_RemoveDura_jist.m` passes `-inCurvature
0.5`, `-inLevel 1.0`, `-inAlignment 0.5`, `-inPressure 0.5`, `-inTarget
0.3`, `-inMask 2.0`, `-inLevel2 0.0`, `-inMaximum 5`, `-inWM3 0.05`, and
`-inGM2 0.02`. The source maps `-inLevel` to alignment force 1.0 and
`-inAlignment` to alpha 0.5; `-inTarget` is the target membership/isovalue
membership 0.3. Do not infer parameter meaning from the flag names.

The dura wrapper topology-corrects WM, then topology-corrects the sum of that
result and GM; the second product is a corrected union, not simply a GM map.
It constructs `inner = 0.5 - wmTC/255` and
`outer = membership - unionTC/255`, runs original JIST `DistanceField` with
caps 15 and 10, then evolves the outer field with
`DynamicPressureTGDM(CONNECT_18_6)`. The P29 Python port covers a reusable
`GenericTGDM`-style evolution kernel. `src/pymacruise/toads_dura.py` now
implements the selected `DynamicPressureTGDM.updatePoint` reference-gradient
alignment penalty and `DynamicPressureCortexProfile.getEffectiveForce` for
one voxel. A committed 5³ fixture calls the selected Java subclass and profile
on interior positive-pressure and mirrored-boundary negative-pressure
examples; Python matches the updates within `5e-10`. This is point-update
evidence only. It does not validate the GenericTGDM evolution loop, topology
veto, narrow-band scheduling, landmine restart, or distance reinitialization.
Although `MaGenericTGDM` and VandyCRUISE `GenericTGDM` share most core update
behavior, they are distinct classes. `setReferenceLevelSet` computes the
reference gradient; `solveBrainMaskSurface` clears the inner and outer bounds,
so the gradient penalty remains active without bound clamping. Full evolution
still needs the selected 4.8/3.6 narrow-band profile, 30 inner updates, five
outer iterations, and signed-distance caps 6/5. See the fixture methodology
and regeneration command in [the dynamic-update reference README](../tests/data/synthetic/dura_dynamic_update/README.md).

The JIST option strings in `MedicAlgorithmToads08` and the matching output
branches in the companion `toads-cruise/src/.../AlgorithmMultipleToads.java`
were checked together; that companion copy uses the lowercase option strings
(`dura removal inputs`, `cruise inputs`) passed by the MATLAB wrapper. It is
the selected adapter reference. The VandyCRUISE source copy uses different
enum strings and is not the right wrapper reference for this entry point.

For TOADS `dura removal inputs`, `AlgorithmMultipleToads` returns five
uncropped products in order: hard classification (`ubyte`), CSF membership,
GM membership, filled-WM membership, and WM-mask membership (four float
volumes). `MedicAlgorithmToads08` multiplies each membership by 255 before
exposing it to the dura-removal wrapper. The MATLAB stage saves all five as
`*_toads1{csf,gm,seg,wmfill,wmmask}.nii.gz` or the corresponding `toads2`
names. Final meaning, encoding, and RAS affines should be verified in the Java
harness instead of inferred from variable names.

### Completed substage: masked robust intensity normalization

`src/pymacruise/toads_preprocessing.py` ports the selected
`MultipleToadPreprocess` initialization of a T1 channel. It applies the
required user mask, reproduces `ImageFunctions.robustMinimum` and
`robustMaximum` using the source 10-bin histogram refinement (fraction 0.01,
two scales), and divides each masked voxel by `Ihigh-Ilow` without subtracting
`Ilow`. This keeps outside-mask zeros at zero, as in the source. The committed
13×11×9 NIfTI phantom includes zero background and bright outliers; all 1,287
normalized float32 values and both bounds match the test-only Java harness
exactly. A non-RAS input test checks the warning, reorientation, and same
normalized RAS values. The implementation rejects a zero robust span with a
clear error rather than propagating division-by-zero NaNs.

Regenerate this substage reference with:

```sh
./.venv/bin/python tests/data/synthetic/toads_preprocessing/make_phantom.py
javac -d build/java-reference tests/java_reference/ToadsIntensityReferenceHarness.java
java -cp build/java-reference ToadsIntensityReferenceHarness > tests/data/synthetic/toads_preprocessing/java_reference.json
```

### Completed substage: exact initial intensity centroids and memberships

`src/pymacruise/toads_initialization.py` ports
`MultipleToadPreprocess.initialCentroids` and the exact
`initialExactMemberships` branch. `AlgorithmMultipleToads` sets
`useApprox=false` in the selected source constructor, so this is the active
initial membership branch for the selected wrapper. Centroids use
`low/(high-low) + atlas_prior`, retaining the source's zero-preserving
normalization convention. Membership scores sum channel-wise squared
distances in Java float order, invert distances above `1e-20f` (otherwise use
`1e20f`), multiply by the exact channel-0 centroid tie count, then normalize.

The committed 6×6×6 two-channel fixture compares all four class centroids and
all 864 output memberships exactly against the standalone Java harness. It
exercises tied centroids, an exact zero-distance sample, and the crop's zero
padding. Because the source multiplies scores by tie-group size after building
the denominator, memberships need not sum to one when centroids tie; the port
preserves that behavior. Atlas descriptor intensity-prior parsing, transformed
atlas memberships, registration, and iterative TOADS segmentation remain open.
Regenerate the reference with:

```sh
javac -d /tmp/pymacruise-toads-init tests/java_reference/ToadsInitializationReferenceHarness.java
java -cp /tmp/pymacruise-toads-init ToadsInitializationReferenceHarness \
  > tests/data/synthetic/toads_initialization/java_reference.json
```

### Completed substage: selected atlas intensity-prior parsing

`load_toads_intensity_priors` in `src/pymacruise/toads_atlas.py` follows the
selected `MultipleToadDeformableAtlas.modalityId` alias mapping and
`Intensity Atlas` descriptor section. It returns float32 priors in declared
class order and in caller-requested channel order. The selected base
T1_MPRAGE profile's 11 values match the standalone Java descriptor harness
exactly. Alias, order, missing-profile, unknown-modality, and malformed-row
cases are covered by `tests/test_toads_atlas_priors.py`. The selected T1
profile is explicit, so this slice reports a missing profile rather than
silently applying the source class's alternate label fallback. It reads only
descriptor text and does not load atlas image arrays.

The source Java atlas reader treats atlas files as headerless byte/float
payloads in descriptor index order. The bundled NIfTIs have been axis-flipped
to RAS and their headers do not necessarily represent the geometry used by
that reader. `load_toads_atlas_data` uses the conversion audit's input axis
codes to reverse only the RAS axis permutation/flips. Its synthetic asymmetric
case checks the restoration without outside source files. A one-off
value-by-value audit against the 12 selected `.raw` files still present in
the supplied source checkout found exact equality for topology and all 11
shape maps (zero mismatches). The raw source files are not bundled; runtime
loading uses only repository NIfTIs and the committed conversion audit.

### Completed substage: original-JIST signed distance field

`src/pymacruise/toads_dura.py::dura_signed_distance_field` adapts the shared
six-neighbor unit-speed fast-marching implementation to the selected dura
stage. It explicitly rounds each intermediate distance update to float32,
matching JIST `DistanceField` storage. The four deterministic fields include
an exact-zero plane, an oblique multi-axis crossing, an interface near a volume
face, and a truncated curved interface. All 648 output voxels match the
original `DistanceField.class` loaded from the supplied JIST-Test JAR exactly.
The comparison does not validate dura mask topology, surface evolution, or the
complete `MedicAlgorithmRemoveDura09` pipeline. See
`tests/data/synthetic/p18_distancefield_jar_reference/README.md` for the exact
JAR/class hashes, fixture recipe, and Java invocation. Focused regression
coverage is in `tests/test_toads_dura.py`.

### Completed substage: Dura dynamic-pressure TGDM evolution

`src/pymacruise/toads_dura.py::dura_tgdm_evolve` now connects the selected
`DynamicPressureTGDM` profile to the reusable `tgdm_evolve` kernel. It applies
the wrapper's gray-matter membership normalization, reference-gradient
alignment, curvature sign, 4.8/3.6/6.0/5.0 distance profile, 30-step outer
passes, time steps, selected 18/6 scheduler, and float32 distance updates.
`tests/test_toads_dura.py::test_dura_tgdm_full_profile_matches_vandy_java_harness`
compares five deterministic outer passes on an 11×10×9 synthetic nested-surface
phantom against the committed Java harness fixture. The five passes schedule
150 steps and produce no landmine vetoes on this case. Signs match at every
voxel and pass. The maximum level-set difference after each reinitialization
is 0.004042, 0.003458, 0.000504, 0.0000117, and 0.00000048 voxel; the final
mean absolute difference is 1.01e-7 voxel. Pre-reinitialization differences
peak at 0.01083 voxel. These are measured discrepancies, not a claim of exact
intermediate parity.

The full-evolution harness compiles the selected VandyCRUISE `GenericTGDM`,
`DynamicPressureTGDM`, and `DynamicPressureCortexProfile` source and runs the
original JIST `DistanceField.class` during iterative evolution. The
reproducible headless setup uses the actual JIST image/distance classes and a
small framework shim for presentation/preferences. Its class SHA-256 is
`92294be2db74351446c3e5898fc35b58497f89ccb38772d4d565885f3f0c35b7`. The
companion `DistanceFieldyh` harness is retained as an auxiliary check; its
full-evolution arrays match the original JIST class exactly on both committed
cases.

The 11×10×9 nested-surface case schedules five 30-step passes with no landmine
restarts. A separate 25×25×25 angularly folded cortex fixture triggers one
source landmine restart in both implementations: Java and Python complete
147 updates. Final output metrics pass the project gates: Dice 0.997923,
symmetric mean surface distance 0.005180 mm, and Hausdorff distance 0.328474
mm. Four voxel signs differ, and intermediate reinitialization fields differ
by up to 0.36 voxel; this is reported explicitly rather than hidden in the
surface tolerances. See
`tests/data/synthetic/dura_tgdm_folded_cortex/README.md` for the fixture and
metrics. Regenerate the actual-JIST references using
`scripts/generate_java_dura_original_jist_reference.sh`.

The selected `run_dura_removal` wrapper now composes the two topology calls,
level-set/distance initialization, dynamic-pressure TGDM, and final output
masking on a checked native T1 RAS grid. A synthetic 9³ composition test runs
one outer TGDM pass and matches its intermediate topology fields to the Java
two-pass harness, then checks assembled products against the validated stage
sequence. Full default five-pass/folded-case wrapper acceptance and integrated
TOADS segmentation remain open.

The reference fixture is `tests/data/synthetic/dura_tgdm_full_evolution/java_reference.json`;
its generator and Java harness are under `scripts/` and `tests/java_reference/`.
It is generated from synthetic fields only. Source SHA-256 values:

- `GenericTGDM.java`: `70f5b21541e8a6c6139f1bd227c4ece12e8818db717f2a78ef50d69fe89f08d1`
- `DynamicPressureTGDM.java`: `53b872a5828539bb51d70db0c048d1fd8162f6b7a6b63176266ef7e2e796b6fc`
- `DynamicPressureCortexProfile.java`: `f85912d621198de297e65d6b9205505b3ac2de72e67e6840f3968c58a850f505`
- `FastMarchingGradient.java`: `cf857c0c62e2c7a219e9e530c863b4e40ae594b55e90a355cd7fdf6e0709d3a3`
- `DistanceFieldyh.java`: `a0df8859367a52a6273c488fed74d3fbdf0bbf5776b6e6cfea15261a5b290811`

The nested full-profile and folded-cortex tests pass. The full suite passed
628 tests with 2 expected RAS warnings in 110.04 seconds on tested commit
`03626d4986c20eb9c3fc934d2b09a7ac75f9a84f`; see
[`docs/reports/full-pytest-2026-09-30-p18-dura-tgdm.md`](reports/full-pytest-2026-09-30-p18-dura-tgdm.md).
The folded final segmentation and surface metrics meet the agreed thresholds;
four voxel signs and the measured intermediate-field differences remain
documented. P18 remains Started because TOADS evolution, registration
optimization, and full end-to-end equivalence remain open.

### Completed substage: Dura wrapper output-mask application

`apply_dura_output_mask` uses NumPy masking for the source wrapper's final
voxelwise output loop. This is generic array work, so the Python port uses
NumPy rather than reimplementing the Java loop. The adapter preserves the
source comparison order (`evolved_level_set >= mask_translation` is outside),
zeros original/WM/GM values there, creates an 8-bit inside mask, and applies
`level_set_translation` only after the threshold decision. The synthetic test
includes values exactly on both sides of the threshold, a translation that
would change the mask if applied first, and verifies inputs are not mutated.
The boundary-order focused test passes. The full suite passed 628 tests with 2
expected RAS warnings in 109.34 seconds on tested commit
`ce00bbcdfe42c69c8136854b310a3d28d267f149`; see
[`docs/reports/full-pytest-2026-09-30-p18-dura-adapters.md`](reports/full-pytest-2026-09-30-p18-dura-adapters.md).
This validates only this voxelwise output adapter; topology correction,
end-to-end input construction, physical headers, and the entire removal
wrapper remain open.

### Completed substage: Dura initial surfaces

`dura_initial_level_sets` forms the inner field as `0.5 - WM/255` and the
outer field as `membership - GM/255`, preserving float32 operation order, then
reuses `dura_signed_distance_field` with the source caps of 15 and 10 voxels.
Its 9×8×7 synthetic fixture runs the exact source formulas and the original
JIST `DistanceField.class` from the supplied JIST-Test JAR. Both membership
level-set arrays and all 1,008 resulting signed-distance voxels match exactly.
The source wrapper SHA-256 is
`1ed180daf4a4334c995ad310828ec91cfd5593413eec17ab15404a421d7805a4`; the
reference JAR `DistanceField.class` SHA-256 is recorded in the fixture. The
generator is `scripts/generate_java_dura_initial_levelsets_reference.sh` and
accepts the original jar through `JIST_TEST_JAR`; the jar is only a development
oracle, not a runtime dependency. The exact initial-field Java comparison
passes. The full suite passed 628 tests with 2 expected RAS warnings in
109.34 seconds on tested commit
`ce00bbcdfe42c69c8136854b310a3d28d267f149`; see
[`docs/reports/full-pytest-2026-09-30-p18-dura-adapters.md`](reports/full-pytest-2026-09-30-p18-dura-adapters.md).
This validates initial arithmetic/distance construction only, not the two
topology-correction algorithms or the complete Dura workflow.

### Completed substage: first Dura WM topology correction

The next implementation slice is only the first topology-correction call in
`MedicAlgorithmRemoveDura09.execute()`: it passes the white-matter membership
field and the original WM mask to `tcWM`, returning `wmTC`. The selected
adapter uses foreground/background connectivities `(6,18)`, object-to-
background propagation, normalized levels 0–1, paint threshold 0.5, and
regularization 0.05. Python should call the existing source-compared
`topology_correct` engine with those explicit settings; its general-purpose
regularization default is different and must not be inherited.

The selected Java classes and Python adapter now match five committed
synthetic fields exactly: spatial-hash distinct ranks; plateaus with values
below 0 and above 1; boundary-touching paint; a thin painted branch; and two
disconnected seeds. All 3,807 output voxels match. Full event tracing on each
case matches 4,648 queue pops and the 1,511 critical-candidate decisions,
including Java's custom priority-tie order. The focused topology and Dura test
modules pass 36 tests. Harness provenance and regeneration are recorded in
`docs/reports/topology-correction-java-comparison-2026-09-29.md`.

This substage does not build `wmTC + raw GM`, call the second `tcGM` correction
(which reuses the original WM mask), initialize/evolve Dura TGDM, or establish
complete Dura or TOADS equivalence.

### Completed substage: Dura WM-plus-GM union and second topology correction

`dura_correct_gray_matter_topology` reuses the generic geometry-checked
`calculate_volume(..., "ADD")` operation to form the wrapper's float32 union
from corrected WM and the original raw GM image. It then runs the second
`tcGM` topology adapter using the same original WM paint mask and the same
Dura-specific 6/18 parameters, returning `gmTC`. The source does not feed
corrected GM into the union or create a separate GM paint mask. Synthetic
tests check float32 addition, the reused mask and settings through the
composition, RAS geometry, and mismatched-grid rejection. Direct Java output
comparison now checks the WM correction, float32 union, and GM correction
against a headless harness that follows the selected wrapper's two-call
sequence; all three arrays match exactly on the committed 9³ fixture. The
harness does not instantiate the full JIST wrapper. TGDM evolution and the
full dura workflow remain open.

### Completed geometry slice: TOADS crop and topology resampling

`src/pymacruise/toads_crop.py` ports the selected source's inclusive crop-bound
calculation, inverse-Rodrigues rotation matrix, zero-padded image crop,
source-style translation update/reversal, and nearest-neighbor categorical
topology sampling. It preserves the Java `[x][y][z]` voxel order and the
source's `dimension / 2.0` center convention. The asymmetric formula-level
fixture exercises mixed axis signs, unequal voxel resolutions, nonzero
rotation, a border-touching image region, atlas labels 2 and 3, and clipped
uncropping. Its test compares bounds and labels exactly, serialized float32
image arrays exactly after parsing, and transform translations within
`2e-7`. The reference is explicitly a standalone transcription of source
formulas; it does not execute the upstream `MultipleToadPreprocess` class.
This slice excludes iterative atlas registration and does not resolve the
separately documented atlas physical-spacing discrepancy. See
`tests/data/synthetic/toads_crop_geometry/README.md` for source hashes, exact
phantom, and regeneration details.

## Porting and validation sequence

1. **Freeze the source path and parameters.** Trace V3 `run_Toads_pipeline.m`,
   `run_Toads.m`, `run_Toads_jist.m`, `run_RemoveDura.m`, and the exact
   `MedicAlgorithmToads08` / `MedicAlgorithmRemoveDura09` implementations.
   Record runtime parameter values, atlas ordering, input/output data types,
   crop and transform conventions, and output products.
2. **Remove preprocessing ambiguity.** Require the caller's brain mask for
   base and plus profiles. Validate it against T1 after RAS normalization;
   neither BET nor SPECTRE is invoked.
3. **Build isolated Java references.** Adapt the test-only harness to accept
   committed synthetic NIfTI fixtures and emit named, checksummed stage
   outputs. Keep Java compilation and reference generation out of the Python
   package runtime.
4. **Port individual TOADS stages.** Translate only algorithm-specific logic
   not already supplied by a compatible Python library. Add a test per
   intermediate stage and preserve source tie-breaking, boundary, and
   convergence behavior.
5. **Port the two-pass orchestration.** Integrate N4, user mask, crop/transfer,
   TOADS pass one, dura removal, and TOADS pass two with explicit RAS geometry
   and restartable stage products.
6. **Validate a larger folded synthetic case.** Exercise narrow tissue bands,
   folded GM/WM interfaces, mask boundaries, and dura-like attachments before
   considering any approved imaging data.
7. **Unblock P34 only after acceptance.** All required outputs must be finite,
   correctly labeled, aligned to the declared RAS reference grid, and meet the
   project Dice/surface tolerances or carry an explicitly accepted deviation.

## Blocker for P34

The Python repository now contains selected TOADS preprocessing, atlas
sampling/registration kernels, optimized CSF/GM/WM per-voxel intensity
likelihood and membership-update terms, selected mixed OBJECT/MASK/OUTLIER
membership passes, and the selected damped centroid update. The selected
Dura removal stages are also composed in a callable RAS-grid adapter and
validated on synthetic inputs. Full TOADS segmentation and accepted
end-to-end Dura equivalence remain incomplete. P34 requires the base workflow
to run and meet its output criteria; selected-method synthetic slices alone
do not unblock it.
Therefore P34 remains **Blocked** until the selected source algorithms are
ported and their reference comparison passes. The P18ae data-term slice is
described in [its Java comparison report](reports/toads-optimized-likelihood-2026-09-28.md),
and P18af validates its PSF-bias input primitive in
[a separate source-formula comparison](reports/toads-psf-ridge-bias-2026-09-28.md).
P18ag validates the six-face smoothing increment in
[another source-formula report](reports/toads-six-face-smoothing-2026-09-28.md).
P18ah validates direct reciprocal-energy normalization in
[its own Java formula report](reports/toads-direct-membership-2026-09-28.md).
P18ai records the direct/neighbor/prior shape-relation gate in
[its Java formula report](reports/toads-shape-relations-2026-09-28.md).
P18aj composes that gate with reciprocal-energy normalization in
[a synthetic source-order comparison](reports/toads-shape-weighted-membership-2026-09-28.md).
P18ak implements and compares one ordered in-place six-face-smoothed sweep in
[a small-grid source-formula report](reports/toads-ordered-membership-2026-09-28.md).
P18as now composes the selected centroid-then-membership sequence for one
optimized-tissue iteration against Java references; see
[the P18as report](reports/toads-optimized-iteration-2026-09-29.md). P18at
matches one topology-aware thinning queue-item transition against the Java
18/6 LUT on nine cases; see [the P18at report](reports/toads-thinning-transition-2026-09-29.md).
P18au ports the selected `BinaryHeap4D` ordering semantics, and P18av composes
that heap with the source thinning-label/boundary-queue initialization. Their
Java comparisons include queue ties, growth, mode changes, full seeded queue
contents, and source array updates; see the [P18au report](reports/toads-binary-heap-2026-09-29.md)
and [P18av report](reports/toads-thinning-queue-2026-09-29.md). Full evolution
still requires dynamic queue updates, repeated thinning and growing, both
wrapper topology-correction stages, deformable atlas registration, convergence,
and dura orchestration. These validated
slices do not replace or establish the full segmentation.

The latest isolated update slices now include optimized-tissue damped
centroids, checked against the actual source `PowerTable` on committed
synthetic data; see the [P18al reference report](reports/toads-optimized-centroids-2026-09-28.md).
The source-compatible optimized-tissue image-variance weights and global
standard deviations are also matched against a Java formula harness in both
field modes; see the [P18am report](reports/toads-optimized-variance-2026-09-28.md).
Assigned-class standard deviations for those optimized tissue types are now
separately checked with the Java formula harness, including source handling
for empty classes; see the [P18an report](reports/toads-optimized-class-variance-2026-09-28.md).
The pairwise centroid-similarity coefficient that consumes those variance
values is also matched for the selected default scale, with source zero
diagonal and IEEE zero-deviation behavior; see the [P18ao report](reports/toads-optimized-similarity-2026-09-28.md).
The voxelwise inverse shape-prior contribution and its half-weighted objective
energy are now matched against the source formula, with the source interior
scan and untouched borders; see the [P18ap report](reports/toads-inverse-shape-prior-2026-09-29.md).
That contribution is now composed into the ordered in-place, six-face
membership sweep and compared end-to-end at formula level on a small synthetic
grid; see the [P18aq report](reports/toads-ordered-prior-sweep-2026-09-29.md).
The optimized CSF/GM/WM per-image terms have now also been connected from
image inputs through the final ordered membership update. Four Java
source-formula cases compare field/PSF configurations, with intermediate and
final outputs meeting the recorded float32 errors; see the [P18ar report](reports/toads-optimized-membership-pass-2026-09-29.md).
P18as also composes one selected optimized-tissue centroid-then-membership
iteration against a Java formula harness; see [the P18as report](reports/toads-optimized-iteration-2026-09-29.md).
P18at validates the local topology-aware thinning transition using the selected
18/6 LUT, P18au validates the Java four-dimensional heap order, and P18av
composes source thinning initialization and boundary-queue seeding. P18aw
validates the post-removal face-neighbor scheduling loop and its source
priority/volume guards. These queue slices preserve two source counter quirks,
recorded in the respective P18av/P18aw reports.

P18ax now composes the pop, source volume/skeleton gates, best-class and
relation target selection, 18/6 topology decision, skeleton/volume update, and
post-removal neighbor scheduling into one queue event. Seven 5³ source-harness
cases match the Java LUT, heap, arrays, float32 class volumes, and remaining
queue exactly. The additional source rule discovered in this composition is
recorded: queue class `k` selects the center membership and volume update, and
the source only gates on `skeleton[x][y][z] == k`, even if the center
segmentation label differs. The [P18ax report](reports/toads-thinning-queue-step-2026-09-29.md)
contains reference hashes and per-case evidence. P18ay now repeats this event
under the source outer `val > -limit` condition and compares every event
snapshot plus the final queue state to Java; see the [P18ay report](reports/toads-thinning-queue-run-2026-09-29.md).
That test begins from a P18ax fixture rather than composing queue initialization
and the repeated loop, and still does not establish full TOADS segmentation.

P18az composes initialization and the repeated event loop on the deterministic
12³ four-class phantom. All 1,000 initial queue entries are processed; the
initial/minimum/final volumes, boundary counts, discrete arrays, and queue
state match the Java harness exactly. A 0.8 volume ratio exercises both the
source 100 floor and ratio-derived minimum branches. See the [P18az report](reports/toads-initialized-thinning-pass-2026-09-29.md).
This validates one source-faithful synthetic thinning pass, while full
registration, other TOADS iteration stages, growing, Dura, and real-data
equivalence remain open.

P18ba ports the separate initGrowingLabels boundary queue: empty-skeleton
availability, six-face class seeds, source float32 priority construction, and
min-heap ordering. Three synthetic cases match exactly with 24 mixed-label
entries, six single-seed entries, and an empty queue. See the [P18ba report](reports/toads-growing-queue-initialization-2026-09-29.md).
P18bb now ports and validates one topology-gated event, including accepted,
critical, and stale queue entries against the selected 18/6 LUT and Java heap;
see [its source-slice report](reports/toads-growing-queue-step-2026-09-29.md).
P18bc composes initialization with queue-to-exhaustion growth on a synthetic
folded-band phantom. Its 448 events and final arrays match the selected Java
harness reference exactly; see [the P18bc report](reports/toads-growing-pass-2026-09-29.md).
This remains source-slice evidence only and does not establish full TOADS
segmentation, atlas alignment, or real-data parity.

### Completed substage: P18bg AlgorithmMultipleToads loop control

The Python `toads_iteration_control` helper matches a Java control-flow
excerpt for the source float32 energy-ratio formula, convergence and maximum
iteration checks, optional growing when stopping, final classification update,
and relation refresh only while continuing. Five cases cover continuation,
convergence, and iteration limit with `doGrowing` enabled and disabled.
Focused test: 1 passed. Full suite: 529 passed with two expected RAS warnings
in 77.99 seconds on implementation commit
`d1e8721e034524062995039724ed93dfbdc18155`; see
[`docs/reports/full-pytest-2026-09-29-toads-loop-control.md`](reports/full-pytest-2026-09-29-toads-loop-control.md)
and
[`tests/data/synthetic/toads_iteration_control/README.md`](../tests/data/synthetic/toads_iteration_control/README.md).
This validates only the loop's post-membership control decisions, not a full
TOADS iteration.

The next TOADS work must continue through source-backed stages and must not
present these formula slices as an end-to-end segmentation.

### Completed substage: P18bh TOADS stopping-iteration classification transfer

The selected outer loop calls `MultipleToadSegmentation.updateClassificationFromLabels()`
when convergence or the iteration cap stops the loop. This source method
replaces each classification byte with `templateLabel[segmentation]` and
returns the count of changed voxels. P18bh isolates that final label transfer;
it does not implement membership-based winner selection, growing, topology,
or the complete outer iteration. The selected source file hash is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`. A
source-excerpt Java harness and the Python mapping match all output bytes and
changed counts on a 3×2×4 phantom; input range/shape checks are covered. The
focused comparison and input checks pass (6 tests). The full synthetic suite
reports 535 passed with 2 expected RAS warnings in 78.80 seconds on
implementation commit `b005798`; see
[`docs/reports/full-pytest-2026-09-29-toads-classification-labels.md`](reports/full-pytest-2026-09-29-toads-classification-labels.md)
and
[`tests/data/synthetic/toads_classification_labels/README.md`](../tests/data/synthetic/toads_classification_labels/README.md).

This validates the source label-transfer method only. The complete
AlgorithmMultipleToads numerical loop and full TOADS evolution remain open.

### Completed substage: P18bi TOADS four-slot closest-label propagation

`toads_closest_labels` in `src/pymacruise/toads_relations.py` ports
`MultipleToadSegmentation.computeClosestLabels`. It seeds slots from the
selected six-face traversal, then propagates missing labels synchronously for
the source's `1 <= t < 5*limit` passes. The Java source-excerpt harness and
vectorized Python function match all 2,016 slot values for three limits on an
asymmetric 9×8×7 six-class phantom. Focused checks: 7 passed. Full suite: 542
passed, 2 expected RAS warnings, in 78.03 seconds on implementation commit
`5f7a1fad86af0be938d8493a17b179b35d6c2e6b`; see
[`docs/reports/full-pytest-2026-09-29-toads-closest-labels.md`](reports/full-pytest-2026-09-29-toads-closest-labels.md)
and
[`tests/data/synthetic/toads_closest_labels/README.md`](../tests/data/synthetic/toads_closest_labels/README.md).
This validates closest-label maps, not relation-map construction or complete
TOADS segmentation.

### Completed substage: P18bj TOADS relation-list construction

`toads_build_relation_list` in `src/pymacruise/toads_relations.py` follows
`MultipleToadSegmentation.buildRelationList`: singleton IDs, first-seen
unordered interface sets, compact interior map IDs, relation membership, and
the class-level union rows. The selected-method Java harness and Python match
all map entries and relation bits for uniform, planar, and folded six-class
synthetic volumes. Focused checks: 7 passed. Full suite: 549 passed with two
expected RAS warnings in 77.03 seconds on implementation commit
`3ebf1248ef891ba5ef772fe4c8232071e2fc412f`; see
[`docs/reports/full-pytest-2026-09-29-toads-relation-list.md`](reports/full-pytest-2026-09-29-toads-relation-list.md)
and
[`tests/data/synthetic/toads_relation_list/README.md`](../tests/data/synthetic/toads_relation_list/README.md).
Only relation-set construction is validated here; diffusion into interior
voxels is a separate step.

### Completed substage: P18bk TOADS relation-map diffusion

`toads_compute_relations` builds the selected relation table and runs the
exact-membership `computeRelations` update. Java and Python match the final
map and pass count across uniform, planar, and folded cases with asymmetric
float32 memberships. The fixtures cover zero limit, zero-support strict ties,
the exclusive iteration cap, and early stop. Focused checks: 5 passed. Full
suite: 554 passed, 2 expected RAS warnings, in 88.89 seconds on implementation
commit `19ae617f0fca849b11a60211758792acedc25ec1`; see
[`docs/reports/full-pytest-2026-09-29-toads-relation-diffusion.md`](reports/full-pytest-2026-09-29-toads-relation-diffusion.md)
and
[`tests/data/synthetic/toads_relation_diffusion/README.md`](../tests/data/synthetic/toads_relation_diffusion/README.md).
This validates exact-membership relation propagation only. Approximate mode,
topology evolution, and full TOADS segmentation remain open.

### Completed substage: P18bl TOADS boundary-distance update

`toads_update_boundary_distances` ports the source's class-wise signed-byte
boundary field update. The selected Java method and Python match every byte on
fresh and prepopulated synthetic arrays. These checks preserve two source
quirks: the diffusion loop clears `change` without restoring it, and the
`min(0, neighbor + 1)` expression cannot expand a zero value. Thus the method
marks both class sides of an interface but does not produce a distance ramp.
Focused checks: 6 passed. Full suite: 560 passed with 2 expected RAS warnings
in 77.53 seconds on implementation commit
`5341266f146ef305c2412f31d846b74dcff01d1d`; see
[`docs/reports/full-pytest-2026-09-29-toads-boundary-distances.md`](reports/full-pytest-2026-09-29-toads-boundary-distances.md)
and
[`tests/data/synthetic/toads_boundary_distances/README.md`](../tests/data/synthetic/toads_boundary_distances/README.md).
This validates boundary-distance state only, not distance-weighted
membership behavior or complete TOADS evolution.

### Completed substage: P18bm source-ordered relation-state refresh

`toads_refresh_relation_state` composes the selected wrapper's
`computeRelations` → `computeBoundaryDistances` → `computeClosestLabels`
sequence. It returns the diffused relation map/table, updated carried
signed-byte boundary state, and four closest-label slots from one segmentation
snapshot, with the selected loops' pass counts. A composed selected-Java
harness and Python match every map/table entry and byte on an asymmetric
folded 9×8×7 phantom with nonuniform float32 memberships and prepopulated
boundary state. Focused checks: 5 passed. Full-suite evidence is recorded in
[`docs/reports/full-pytest-2026-09-29-toads-relation-state-refresh.md`](reports/full-pytest-2026-09-29-toads-relation-state-refresh.md)
for 565 passed and 2 expected RAS warnings in 77.82s on implementation commit
`ae7e5558b301c22f1ebe4d334a1fda7e90b6a711`. This validates selected state
refresh only; approximate-membership mode and TOADS segmentation evolution
remain open.

### Completed substage: P18bq TOADS ordinary OBJECT intensity energy

`toads_object_intensity_energy` reproduces the selected ordinary OBJECT
branch of `computeExactMemberships`: source-order modality-weight ×
class-variance × image-variance scaling and squared intensity residuals, with
optional field correction and PSF bias. NumPy vectorizes across the synthetic
volume while each modality/class operation retains float32 behavior. Selected
Java and Python match every output bit for raw and corrected three-modality,
four-class cases, with and without a PSF-bias array. Focused checks: 2 passed.
Full suite: 578 passed, 2 expected RAS warnings, in 78.57s on implementation
commit `dd26e5ce978514c0df75f022bf8a7bf370d1ea71`; see
[`docs/reports/full-pytest-2026-09-29-toads-object-intensity-energy.md`](reports/full-pytest-2026-09-29-toads-object-intensity-energy.md)
and [`tests/data/synthetic/toads_object_intensity_energy/README.md`](../tests/data/synthetic/toads_object_intensity_energy/README.md).
This is one source data term, not a full class-energy or membership update.

### Completed substage: P18br ordinary OBJECT membership pass composition

`toads_object_membership_pass` composes the P18bq data term with the existing
inverse-shape prior and `toads_ordered_membership_sweep`. It reuses the
source-checked Python helpers for prior terms, in-place six-face smoothing,
relation weighting, reciprocal conversion, normalization, and objective
accumulation. The selected Java composition harness and Python match every
membership bit and objective float on a 4×4×4 synthetic volume with two
modalities and three classes, in both raw and field-corrected/PSF-bias modes.
Focused checks: 2 passed. Full suite: 569 passed, 2 expected RAS warnings, in
77.32s on implementation commit `e01b03e838fa9330f7487c9feb31260a3e1883dd`;
see [`docs/reports/full-pytest-2026-09-29-toads-object-membership-pass.md`](reports/full-pytest-2026-09-29-toads-object-membership-pass.md)
and [`tests/data/synthetic/toads_object_membership_pass/README.md`](../tests/data/synthetic/toads_object_membership_pass/README.md).
This is the ordinary OBJECT pass only; lesion/outlier mixtures, topology
gates, and outer TOADS evolution remain open.

### Completed substage: P18bs TOADS MASK intensity energy

`toads_mask_intensity_energy` implements the selected MASK likelihood using
the source `maskDistance(image, 0)`: exact positive and negative zero map to
zero cost, and all finite nonzero intensities map to unit cost. The remaining
per-modality and per-class weights accumulate in Java float32 order. A
selected Java formula harness and Python match every energy bit on synthetic
multimodal inputs. Focused checks: 2 passed. Full suite: 571 passed, 2
expected RAS warnings, in 78.37s on implementation commit
`8cc0d1ecc491a187b114b613930c6443d7f94953`; see
[`docs/reports/full-pytest-2026-09-29-toads-mask-intensity-energy.md`](reports/full-pytest-2026-09-29-toads-mask-intensity-energy.md)
and [`tests/data/synthetic/toads_mask_intensity_energy/README.md`](../tests/data/synthetic/toads_mask_intensity_energy/README.md).
This is the MASK data term only; mixed object types and full TOADS evolution
remain open.

### Completed substages: P18by, P18bz, and P18ca exact-membership compositions

`toads_object_mask_outlier_membership_sweep` composes class and scalar
OUTLIER energies, the source-ordered smoothing and relation factors, inverse
reciprocals, normalization of both class and outlier maps, and objective
accumulation. P18bz adds and validates the selected atlas inverse-prior and
half-weighted objective terms. P18ca adds
`toads_object_mask_outlier_membership_pass`, which dispatches the selected
OBJECT residual, MASK mismatch, and paired OUTLIER intensity terms into that
sweep. The selected Java formula harness matches every tested float32 bit on
the committed synthetic phantoms, including field correction, PSF bias, and
atlas priors. Full-suite evidence and exact implementation commits are linked
from P18by–P18ca in [`plan.md`](../plan.md). These compositions still do not
perform TOADS deformation, topology-gated growth/thinning in the complete
outer loop, or final segmentation transfer.

### Completed substages: P18cy optimized-tissue dispatch and P18cz Rayleigh distance

P18cy extends fixed-input `computeExactMemberships()` composition to
interleaved OBJECT, MASK, OPTIMIZED_CSF, OPTIMIZED_GM, and OPTIMIZED_WM
classes. The pinned Java harness covers category-specific PowerTable
exponents, field and PSF modes, direct/neighbor/prior relations, smoothing,
atlas priors, objective order, and unchanged borders. Focused checks: 7
passed. P18cz implements only `rayleighDistance(float val, float sig)`;
17 Java/Python cases match float32 bits, including Java-compatible
nonpositive returns and IEEE edge outcomes. Focused checks: 2 passed. The
full suite reports 736 passed with 2 expected RAS warnings in 168.07 seconds
on implementation commit `a203b7ec912046b07017de34b63cb595a5f887af`; see
[`docs/reports/full-pytest-2026-09-30-toads-rayleigh-lesion-energy.md`](reports/full-pytest-2026-09-30-toads-rayleigh-lesion-energy.md)
and the P18cy/P18cz fixture provenance in `tests/data/synthetic/`.
Optimized dispatch still excludes LESION, OUTLIER, OUTMASK, and Rayleigh
class mixing. The Rayleigh primitive does not estimate centroids or route a
class family. Neither slice closes the TOADS iteration or workflow.

### Completed substage: primary TOADS output image creation (P18cyy)

`toads_output_io.py` serializes the selected primary `AlgorithmMultipleToads`
result modes (`hard_segmentation`, `full_segmentation`, `Cruise_input`,
`DuraRemoval_input`, and `thalamus_seg`) from already-restored RAS arrays.
It preserves the source `_Toads_seg`, `_Toads_mems`, `_csf`, `_gm`, `_wmfill`,
`_wmmask`, `_thalamus`, and optional `_Toads_field(s)` suffixes, with uint8
labels and float32 membership/input/field data. A shared image helper applies
the RAS T1 geometry; files are saved as `.nii.gz`. The debug `all_images`
branch and crop restoration are not included. The focused output/adapter suite
passed 33 tests in 1.47 seconds on commit
`86ba25a79d889508c2f710c01ba1bf9535b72032`. This validates serialization only,
not full TOADS workflow equivalence.

### Completed substage: Dura-output wrapper-buffer to NIfTI bridge (P18cyz)

`save_toads_dura_wrapper_outputs` consumes the selected wrapper's flattened
classification and four 255-scaled float buffers. It reshapes the Java
x-fastest buffer order using Fortran indexing into the NIfTI `[x,y,z]` grid,
preserves UBYTE labels and float32 outputs, then uses the P18cyy RAS NIfTI
writer and selected suffixes. The asymmetric 2×3×4 phantom and committed Java Dura-output fixture check
every voxel after compressed save/reload. The Java fixture also passes through
wrapper scaling and Fortran-order reconstruction. Focused six-file
output/adapter suite: 40 passed in 0.57 seconds on implementation commit
`60ad9586cf8341f9312ed6fa22f4c1689f943e23`. This bridges the wrapper's final
Dura products only; it does not demonstrate a complete TOADS workflow.

### Completed substage: one-channel TOADS alignment startup composition (P18cyaa)

`src/pymacruise/toads_alignment_startup.py` connects the selected
`AlgorithmMultipleToads` setup calls through creation of the initial exact
membership volume: mask and robust-span normalization, crop-bound discovery,
image cropping with the source padding and transform update, topology sampling
onto the crop grid, intensity-prior centroid initialization, and the exact
membership initializer. The module returns each intermediate and an explicit
source-stage order so registration can consume the result without guessing
about geometry or shape.

The adapter is currently the selected one-channel base-TOADS path. It accepts
atlas topology already expressed in RAS voxel-axis order and takes its spacing
from the source descriptor. The focused startup/preprocessing/crop/initialization
group passes 23 tests in 0.64 seconds on implementation commit
`02ac7639900c4feaa58e9d15fbeb2482c5d44daf`. The returned state separately
records the original T1 orientation codes while its image geometry describes
the normalized RAS grid. The component functions have their
own selected-source evidence; this composition does not have an independent
Java end-to-end reference and ends immediately before rigid atlas registration.
It does not establish registration, iterative segmentation, or full TOADS
equivalence.

### Completed substage: initial exact-membership atlas-pyramid call (P18cyab)

`register_toads_alignment_startup` continues P18cyaa by passing its exact
membership volume, crop-adjusted transform, native subject spacing, and
descriptor template spacing to `register_toads_atlas_rigid_pyramid`. It
excludes class zero (the selected MASK/background class) from the optimizer,
as `MultipleToadDeformableAtlas` does. Shape maps are supplied in descriptor
class order; atlas loading and the preceding `updateShapeRegistration` state
setup are not folded into this adapter.

The focused startup and registration-sampling group passes 45 tests in 1.12
seconds on implementation commit `23144f010e55fbf9509fb70f108510ecf8742771`.
This checks composition with an asymmetric synthetic atlas and exercises the
registration kernels that have separate selected-Java evidence. It is not an
independent Java comparison of the complete startup-to-registration chain and
does not validate TOADS iteration or full workflow equivalence.


### Completed substage: RAS-oriented atlas-array loading (P18cyac)

`load_toads_atlas_ras_data` now provides algorithm-facing NIfTI arrays in
canonical RAS voxel-axis order, with warnings and source-axis metadata when an
asset needs reorientation. It checks finite 3D data and a shared atlas grid,
requires integer topology labels, and returns the descriptor template spacing.
`toads_atlas_shapes_in_structure_order` stacks named shape maps as
`(x,y,z,classes)` in the descriptor structure order. The earlier
`load_toads_atlas_data` remains the source-index-order reader for Java
comparisons; these APIs have different, documented purposes.

The focused atlas/startup/registration group passed 32 tests in 1.49 seconds on
implementation commit `a1fa361`. Its synthetic asymmetric LPS fixture verifies
three warnings, exact flip behavior, retained descriptor spacing, orientation
metadata, and class order; see [the P18cyac report](reports/toads-atlas-ras-loader-2026-10-02.md).
This loader is not yet assembled with descriptor intensity-prior lookup and the
startup/registration adapters, and it does not validate atlas registration,
TOADS iteration, or real-image behavior.


### Completed substage: atlas-assets to initial-registration composition (P18cyad)

`register_toads_from_atlas_assets` accepts a loaded TOADS atlas or descriptor
path, loads its RAS topology/shapes and requested explicit modality prior,
prepares the one-channel T1/mask startup, orders the atlas maps according to
the descriptor, and invokes the selected exact-membership rigid pyramid. The
levels and maximum iterations are explicit inputs, so this helper does not
quietly invent wrapper configuration.

A compact synthetic two-class descriptor and NIfTI atlas/T1/mask exercise that
path through one scale-2 optimizer level. The focused group passed 13 tests in
1.24 seconds on implementation commit
`e80ef81a40195b195eae32c31a4033b8e95262c3`; see [the P18cyad report](reports/toads-atlas-startup-composition-2026-10-02.md). This composition does not
recreate all internal transform-model state transitions of Java
`initShapeRegistration`/`updateShapeRegistration`, compare end-to-end registration
against Java, or run TOADS iterations. The next integration slice should pin
those selected wrapper state transitions before increasing the claim.


### Completed substage: selected rigid-registration setup handoff (P18cyae)

`register_toads_alignment_startup` now exposes a `ToadsRigidRegistrationSetup`
record for the selected `AlgorithmMultipleToads` source sequence. It preserves
the exact-membership array reference, null `best` candidates, `approx=classes`
with `useApprox=false`, six-element zero initialization followed by the
preprocessing transform, rigid mode, iteration/level inputs, subsample factor
2, and image center. The existing numerical pyramid consumes this recorded
state's memberships, transform, level count, and iteration count.

Focused startup/registration tests passed 27 cases in 1.48 seconds on
implementation commit `845eecfa43b9862ce81e7d9322aeedccf4593b7b`. The report
pins the selected source files and lists the values tested. This is a setup-state
translation plus synthetic wiring check, not execution of the Java wrapper or
full TOADS comparison; subsequent algorithm steps remain open.


### Bundled atlas payload smoke check

The selected `cruise-atlas-10obj-toads2008.txt` and all its compressed maps were
loaded with the Python RAS data API at implementation commit
`845eecfa43b9862ce81e7d9322aeedccf4593b7b`. The 11 class maps stacked in
descriptor order to 172×216×164×11; topology labels were exactly the descriptor
set, the T1 prior vector had 11 finite values, and the descriptor spacing was
0.9375 mm. The operation materialized 542,270,592 bytes of arrays. It read
bundled atlas resources only, not subject MRI, and it is not a segmentation or
registration validation. Details and the command are in [the smoke report](reports/toads-bundled-atlas-ras-smoke-2026-10-02.md).


### Corrected signed RAS spacing handoff (P18cyaf)

The startup-to-registration adapter now passes
`ToadsRASGeometry.signed_image_resolutions_mm` unchanged. The selected
`MultipleToadPreprocess.getSignedImageResolutions()` values for RAS voxel axes
are `(-sx, -sy, +sz)`; the previous `abs(...)` conversion discarded source
coordinate direction in the registration map. Descriptor template spacing
continues to be passed separately.

The focused startup and registration group passed 43 tests in 1.25 seconds on
implementation commit `0d15a143d4f8321266e82e8272c5910bac92fb56`. It includes
the Java-reference signed-landmark mapping test and a new adapter-level
assertion that the exact signed tuple reaches the setup record. The source class
hash and scope are recorded in [the P18cyaf report](reports/toads-signed-ras-spacing-2026-10-02.md).


### RAS geometry recorded at registration setup (P18cyag)

The inspectable registration setup also records the `setImageInfo` contract:
subject and template dimensions, selected MIPAV orientation enums `(0, (2,3,5))`,
RAS direction signs `(-1,-1,+1)`, signed subject spacing, descriptor template
spacing, and both voxel-grid centers. These values are derived from the RAS
startup geometry and the actual membership/template-prior arrays consumed by
the optimizer.

The folded asymmetric synthetic startup uses subject spacing `(0.8,1.2,1.5)` mm
and 1-mm atlas spacing. The focused startup and registration tests pass 43 cases
in 1.20 seconds on tested tree commit
`cbd0887d248502e07cc04668d39d2206b3b38527` (geometry implementation commit
`f6b8848bea1a3d76e9ea1d2c3e1406bc28b5bce8`). Existing Java-reference
signed-coordinate landmark tests are included; the new assertions specifically
check the setup geometry. This does not validate later registration levels or
full TOADS behavior. See [the P18cyag report](reports/toads-registration-ras-geometry-2026-10-02.md).


### Selected wrapper defaults for atlas startup (P18cyah)

`register_toads_from_atlas_assets` defaults to the selected wrapper profile:
`T1_MPRAGE`, four registration pyramid levels, and 50 iterations. Callers can
continue to specify smaller settings for fast synthetic slices; these are
explicit overrides rather than silent changes to the selected profile. The
source call and focused tests are recorded in [the P18cyah report](reports/toads-wrapper-defaults-2026-10-02.md).


### Atlas-backed initial alignment through template refresh (P18cyai)

`register_toads_from_atlas_assets` now continues the selected source route after
rigid atlas fitting: it passes the registered transform to the existing
post-registration refresh, restores the normalized T1 to its native voxel grid,
recomputes crop bounds with the atlas background label, crops under the fitted
transform, and transforms topology onto the complete recropped grid. The mask
label comes from the unique `mask` structure in the atlas descriptor; the
background ratio is a caller-visible option with source-default zero.

The compact synthetic atlas composition plus its focused startup, refresh, and
registration tests pass 28 cases in 1.87 seconds on implementation commit
`4adb642a69011e15a8b91476ae966ac68e3ade9e`. The crop/template refresh component
has separate selected-Java reference checks, but this new descriptor-to-refresh
chain is only synthetic wiring evidence. It ends before segmentation object
construction, reinitialization memberships, single-scale alignment, and TOADS
iteration. See [the P18cyai report](reports/toads-registration-template-refresh-2026-10-02.md).


### Atlas-backed post-refresh exact registration (P18cyaj)

The atlas-assets startup now continues after the post-registration crop and
topology refresh. It initializes exact memberships on the refreshed crop,
recomputes transformed shape bounds, excludes the atlas background class, and
performs the selected scale-one rigid `registerShapes()` update. Subject
spacing remains the signed RAS tuple and template spacing remains the
descriptor value. The selected source default is five iterations; the compact
synthetic integration uses an explicit one-iteration override.

The focused group passed 44 tests in 6.37 seconds on implementation commit
`ee350877e61f97adb050d185bfb34f7bfd067ade`, and the full synthetic suite passed
1,049 tests with 11 expected RAS warnings in 1,682.82 seconds on that same
commit. See the [P18cyaj report](reports/toads-post-refresh-registration-2026-10-02.md)
and [full-suite report](reports/full-pytest-2026-10-02-toads-atlas-registration.md).
This validates selected synthetic composition and the separately checked
registration kernel, not a Java end-to-end registration comparison or TOADS
segmentation iteration. The next driver slice is segmentation construction
and initialization of centroids, intensity ranges, lesion centroid, and the
wrapper-selected atlas-aware exact/approximate membership overload.


### Segmentation initialization state (P18cyak)

`register_toads_from_atlas_assets` now follows the selected post-registration
initialization calls. It copies initialized centroids into the segmentation
centroid and prior-centroid state, records the source robust maximum and
intensity scale, maps the descriptor lesion prior into normalized intensity
space, and computes a fresh exact membership array on the refreshed image
crop. The base wrapper explicitly uses `useApprox=false` and sets `approx` to
the atlas class count, so its atlas-aware membership overload selects the
exact path. If the descriptor has no `Lesions` section or does not list a
requested modality, `getLesionPriors` supplies zero; the port follows that
fallback.

The focused atlas/startup/initialization group passed 59 tests in 6.93 seconds
on implementation/test commit `04e7e9b346cbca0ae7285e4cf2735f6e3c56d7b7`. See the
[P18cyak report](reports/toads-segmentation-initialization-2026-10-02.md).
This composes previously checked source operations and does not validate the
TOADS constructor, variance, correction-field, or iterative membership stages.


### Transformed atlas through constructor and relation state (P18cyal)

The source `MultipleToadSegmentation` constructor maps the transformed template
labels to descriptor-order class indices and seeds segmentation, skeleton,
relation, availability, volume, and boundary-distance arrays. The new adapter
recreates that class-index state on the refreshed RAS crop, then supplies it
and the freshly initialized exact memberships to the already source-checked
relation/boundary/closest-label sequence with the selected `limit=3`.

The focused synthetic atlas/startup and relation group passed 32 tests in 9.42
seconds on implementation commit
`cd8d9bdf1854e2c0e79118780071c0868de7c8e1`. It checks descriptor order,
unknown labels, constructor array values and shapes, and the relation-state
outputs. This is wiring evidence over synthetic data, not an independent Java
comparison of the whole segmentation constructor. Variance and initial
membership-energy computation remain outside this slice. P18cyam now composes
initial variance state; the selected unity correction-field and first
membership/energy stages remain open.


### Initial exact class variances (P18cyam)

After membership seeding and relation-state initialization, the selected
wrapper calls `segmentation.computeVariances()`. The startup adapter uses the
already source-checked assigned-class standard-deviation kernel for the
selected MASK and optimized CSF/GM/WM classes. This branch scans the complete
image grid in source x/y/z order and uses membership squared plus `1e-30f` as
its weight. The selected correction fields begin as unity, so the adapter
avoids allocating and multiplying a redundant unity field.

The focused atlas/startup, relation, and Java-reference variance group passed
38 tests in 7.78 seconds on implementation commit
`1aaa5419dc4a97086f6f80cd3e1ceb18603b58a2`; see the [P18cyam
report](reports/toads-initial-variance-2026-10-02.md). This is a component
composition, not the initial membership/energy pass or a full TOADS run.


### Initial image-wide variance weights (P18cyan)

The selected `computeMemberships()` startup sequence next computes the
image-level variance weights used by the optimized tissue membership terms.
The adapter carries the normalized image, seeded exact memberships, initialized
centroids, and class deviations into the already source-checked image-variance
kernel. It excludes MASK, processes the full grid, and preserves source
float32 accumulation and zero-residual behavior. Initial correction fields are
unity, so no extra field array is introduced.

The focused synthetic atlas/startup and Java-reference kernel group passed 45
tests in 8.01 seconds on implementation commit
`5a4aca2567bcfa138020ab3d410337e1b55cee20`; see the [P18cyan
report](reports/toads-initial-image-variance-2026-10-02.md). This slice does
not yet execute the initial exact membership/energy sweep.


### Initial centroid-similarity prelude (P18cyap)

The exact base-TOADS membership routine calls image-variance calculation, then
computes pairwise centroid similarities inline before its voxel loop. The
startup adapter carries the selected non-lesion descriptor order, initial
centroids, class deviations, and normalized image-variance weights through the
already Java-checked float32 similarity kernel. It preserves the descriptor
class matrix including MASK, default `priorScale=0.05`, and zero diagonal.
Lesion-specific pair dispatch and voxelwise membership/energy remain separate.

The focused atlas-startup, optimized-similarity Java-reference, and pair
similarity group passed 8 tests in 7.82 seconds on implementation commit
`8663b867a912c5a5ff7a8ac0b009b75512e0fd15`; see the [P18cyap
report](reports/toads-initial-similarity-2026-10-02.md). The composition check
uses a two-class synthetic atlas and does not yet validate an 11-class complete
startup or execute the exact membership pass.


### Initial statistics through exact membership pass (P18cyaq)

The selected startup now has a bounded composition that computes assigned-class
deviations, image-wide variance weights, and centroid similarities in source
order, then passes them to the exact membership/energy sweep. Atlas shape
sampling, relation maps, closest labels, exponents, modality weights, and
source priors remain explicit fixed inputs at this seam. This composition
checks the stage wiring; it does not make those upstream maps part of the
registration startup or claim a combined Java result.

The focused synthetic RAS and selected-kernel reference group passed 18 tests
in 4.40 seconds on implementation commit
`6fd8b12b6b048f3cf8f25bcb041b0d073542e13b`; see the [P18cyaq
report](reports/toads-initial-statistics-membership-2026-10-02.md). The exact
membership kernel and component statistics have separate Java-reference
checks; atlas transform sampling and iterative TOADS remain open.


### Initial statistics, exact membership, and first thinning (P18cyar)

The selected startup composition now takes the initial RAS membership and
segmentation state through relation/closest-label refresh, class deviations,
image variance, centroid similarity, exact membership/energy, and the first
thinning call. The atlas shape samples and registration state are still fixed
upstream inputs. The focused adapter test compares intermediate relations and
final thinning arrays to direct calls of the same components; it does not claim
a combined Java wrapper result.

The focused synthetic group passed 14 tests in 2.87 seconds on implementation
commit `ada020ac149e2889ede1195202101bc702adcc28`; see the [P18cyar
report](reports/toads-initial-statistics-thinning-2026-10-02.md). Selected Java
references remain at the component level, and topology-changing behavior on a
large volume remains unvalidated in this composition.


### Transformed atlas shape volume on the refreshed RAS crop (P18cyas)

The atlas startup adapter now turns descriptor-order template shape maps into
a class-last volume on the post-registration subject crop. It uses the
selected final transform and per-class current-transform bounds, signed RAS
spacing, source interpolation-term order, shape-slope recovery, and the
class-zero outside-bounds fallback. The existing source point sampler remains
the reference for individual values. The new volume is available in the
startup result, while the exact membership composition still accepts shape
volumes as explicit inputs.

The focused atlas-startup, shape-sampler, recrop, and registration group
passed 35 tests in 7.77 seconds on implementation commit
`8088396d727e38187bc8bddfdd6db08d01198dda`; see the [P18cyas
report](reports/toads-shape-volume-2026-10-02.md). The volume comparison uses
`3e-7` absolute tolerance against the pointwise adapter and is not an
independent Java whole-volume reference.


### Registered startup through refreshed statistics and exact membership (P18cyat)

`compute_registered_toads_initial_membership` carries the recropped image,
initialized memberships and centroids, constructor segmentation, descriptor
class order, closest-label map, and post-registration transformed shape
samples from one atlas-startup result into the selected statistics and exact
membership composition. The channel-first RAS image is moved to the
class-last convention used by the numerical kernels. Exponents, tissue
modality weights, atlas shape/smoothing weights, and effective coefficients
remain explicit inputs because the current startup result does not resolve
those settings. This adapter does not estimate the correction field, run
TOADS iteration, or establish a whole-wrapper Java comparison.

The focused synthetic group builds an actual selected 11-class RAS atlas,
runs atlas registration, and compares the resulting adapter output with a
direct call to the exact-membership kernel. Six tests passed in 17.79 seconds;
the full suite passed 1,060 tests with 11 expected RAS warnings in
1,862.37 seconds on implementation commit
`dcb707fec0425c3e28fe7b20bf2b2c5822a061be`. See the [P18cyat
report](reports/toads-registered-initial-membership-2026-10-02.md) and [full
suite report](reports/full-pytest-2026-10-02-registered-startup-membership.md).
This establishes synthetic integration only, not a combined Java wrapper
oracle or complete TOADS iteration.
