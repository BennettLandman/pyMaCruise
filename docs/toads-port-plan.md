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
0.9375 mm. Its source MIPAV XML and the bundled `.nii.gz` instead declare
1.0 mm. The selected companion `toads-cruise/src/.../MultipleToadDeformableAtlas.java`
parser reads dimensions and resolutions from the text descriptor, then its
loader reads exactly `Nx*Ny*Nz` bytes from the referenced file and does not
use that image header's affine. The port
therefore retains both sources in `load_toads_atlas`: descriptor spacing is
the algorithm's declared atlas-grid spacing, NIfTI spacing is reported as an
explicit warning, and no atlas payload or metadata is silently rescaled.
The root descriptor's filenames were changed to `.nii.gz` for the Python
repository's storage contract; the Java source reader expects raw payloads.
The Python source-index array loader is implemented and exact against the 12
available raw source maps. A future test-only Java comparison adapter must
serialize these arrays in the source primitive type and byte order. Physical
atlas geometry is still an open P18 validation item.

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
