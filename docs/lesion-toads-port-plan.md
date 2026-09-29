# LesionTOADS and MaCRUISEplus port plan

**Status:** P39 is Started. P39a (numeric intensity initialization), P39b
(atlas-prior parsing), P39c (capped six-neighbor distance support), P39d
(topology NIfTI loading), and P39e (shape-prior NIfTI loading) are validated
on synthetic/reference data. The atlas-loader rows validate bundled
compressed-NIfTI handling only, not the Java raw-byte loader, registration,
transformed-prior sampling, or lesion algorithm. This document tracks
algorithm scope and source contracts; it does not claim a full LesionTOADS
port.

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

### Started substage: P39g transformed shape-prior derivatives

LesionTOADS registration differentiates each sampled prior in template
coordinates with `ImageFunctions.linearInterpolationXderivative`,
`Yderivative`, and `Zderivative`, then composes those values with the rigid
coordinate Jacobian. The selected helper already has a source-comparison
fixture in the TOADS registration tests; P39g makes that derivative available
through the LesionTOADS sampler while preserving its own prior-value and
input-bound behavior. This does not implement registration updates or
segmentation.
