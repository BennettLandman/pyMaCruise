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

### Started substage: Dura dynamic-pressure TGDM evolution

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

Attached dura, outer masks, the two pre-evolution topology corrections, and
post-evolution end-to-end wrapper behavior remain open before P18 can be Done.

The reference fixture is `tests/data/synthetic/dura_tgdm_full_evolution/java_reference.json`;
its generator and Java harness are under `scripts/` and `tests/java_reference/`.
It is generated from synthetic fields only. Source SHA-256 values:

- `GenericTGDM.java`: `70f5b21541e8a6c6139f1bd227c4ece12e8818db717f2a78ef50d69fe89f08d1`
- `DynamicPressureTGDM.java`: `53b872a5828539bb51d70db0c048d1fd8162f6b7a6b63176266ef7e2e796b6fc`
- `DynamicPressureCortexProfile.java`: `f85912d621198de297e65d6b9205505b3ac2de72e67e6840f3968c58a850f505`
- `FastMarchingGradient.java`: `cf857c0c62e2c7a219e9e530c863b4e40ae594b55e90a355cd7fdf6e0709d3a3`
- `DistanceFieldyh.java`: `a0df8859367a52a6273c488fed74d3fbdf0bbf5776b6e6cfea15261a5b290811`

At the time this first harness slice was added, the focused Dura test module
passed 11 tests. The later folded-cortex case is separately gated against the
owner-approved Dice and mesh-distance thresholds. P18 remains Started pending
the whole-evolution original-JIST comparison and remaining wrapper stages.

### Started substage: Dura wrapper output-mask application

`apply_dura_output_mask` uses NumPy masking for the source wrapper's final
voxelwise output loop. This is generic array work, so the Python port uses
NumPy rather than reimplementing the Java loop. The adapter preserves the
source comparison order (`evolved_level_set >= mask_translation` is outside),
zeros original/WM/GM values there, creates an 8-bit inside mask, and applies
`level_set_translation` only after the threshold decision. The synthetic test
includes values exactly on both sides of the threshold, a translation that
would change the mask if applied first, and verifies inputs are not mutated.
The focused Dura test module then passed 14 tests. This validates only this
voxelwise output adapter; topology correction, end-to-end input construction,
physical headers, and the entire removal wrapper remain open.

### Started substage: Dura initial surfaces

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
oracle, not a runtime dependency. The focused Dura test module now passes 16
tests, including the nested and larger folded-cortex TGDM comparisons.
This validates initial arithmetic/distance construction only, not the two
topology-correction algorithms or the complete dura workflow.

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

The Python repository currently contains TOADS-like adapter tests but has no
TOADS or dura-removal implementation in `src/pymacruise/`. P34 requires a real
base workflow and cannot be validated by wiring synthetic TOADS-like inputs
through substitute behavior. Therefore P34 remains **Blocked** until the
selected source algorithms are ported and their reference comparison passes.
