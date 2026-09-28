# TOADS, SPECTRE, and BET implementation plan

**Status:** P18 is in progress. P34 base orchestration is blocked until the
TOADS route below is implemented and validated. The user-supplied brain mask
decision removes BET/SPECTRE skull stripping from the runtime contract; it does
not establish algorithmic equivalence to either tool.

## Selected implementation decisions

| Component | Decision | Runtime requirement | Validation evidence still required |
|---|---|---|---|
| TOADS tissue segmentation | Port the selected VandyCRUISE `MedicAlgorithmToads08` behavior and the MaCRUISE two-pass adapter contract into Python. Reuse NumPy/SciPy for array and image operations where they match source semantics; keep atlas alignment, iteration order, labels, masks, stopping rules, and wrapper settings source-derived. Do not substitute a generic classifier. | No MIPAV/JIST Java runtime in the finished package. Bundle all required atlas/configuration inputs under `data/atlases/` as regular Git files; image payloads stay `.nii.gz`. | Build a Java comparison harness around the supplied implementation. First compare intermediate/terminal fields on deterministic synthetic intensity, class-boundary, bias, and mask cases. Then compare final GM/WM-fill/WM-mask labels and geometry. Record voxel Dice per label, exact label agreement where deterministic, and all numeric deviations. |
| Dura removal between TOADS passes | Port the active `MedicAlgorithmRemoveDura09` stage and reproduce its selected placement between the two TOADS runs. Do not replace it with generic hole filling or morphology. | No separate external executable. | Use synthetic attached dura, disconnected outside tissue, thin bridges, and no-dura cases. Compare Java outputs and preserve exact masks/topology where the source is deterministic. |
| BET skull stripping | Exclude BET from the selected Python workflow. The user supplies the brain mask, which must be binary, nonempty, and on the T1 physical grid. | No FSL BET executable or Python BET wrapper. | Test binary/value validation, RAS reorientation warning, and exact shape/affine agreement against T1. This is an input-contract substitution, not a claim that the provided mask reproduces BET. |
| SPECTRE skull-strip cleanup | Exclude SPECTRE from the selected runtime path because the caller supplies the required brain mask. Retain its source only as historical behavior evidence. Do not install or require FreeSurfer. | No SPECTRE plugin runtime and no FreeSurfer. | Document the workflow difference. Compare the port with MaCRUISE references supplied with masks; do not claim SPECTRE parity. If future evidence shows another SPECTRE product is consumed beyond the brain mask, add it as a separately approved scope item. |

## Source inventory and available material

The active TOADS implementation is present in the supplied companion
`toads-cruise` tree, rather than the MIPAV checkout. The key code reviewed so
far is:

- `toads-cruise/src/edu/jhu/ece/iacl/plugins/classification/MedicAlgorithmToads08.java`:
  MIPAV/JIST wrapper, selected output type, isotropic sinc resampling, and
  wrapper parameter defaults/forwarding.
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
1.0 mm. The Java `MultipleToadAtlas` parser reads dimensions and resolutions
from the text descriptor, then its loader reads exactly `Nx*Ny*Nz` bytes from
the referenced file and does not use that image header's affine. The port
therefore retains both sources in `load_toads_atlas`: descriptor spacing is
the algorithm's declared atlas-grid spacing, NIfTI spacing is reported as an
explicit warning, and no atlas payload or metadata is silently rescaled.
The root descriptor's filenames were changed to `.nii.gz` for the Python
repository's storage contract; the Java source reader expects raw payloads,
so a future test-only Java comparison adapter must decompress and serialize
the exact voxel arrays in the source byte order. Physical atlas geometry is
still an open P18 validation item.

## Selected MaCRUISE wrapper contract

The selected `run_Toads_pipeline.m` sequence is N4 correction, skull-mask
preparation, TOADS pass one (`dura removal inputs` output), dura removal, and
TOADS pass two (`dura removal inputs` is the default final mode). The user
brain mask replaces skull stripping; it must mask the N4 image supplied to
TOADS. We still need to verify whether the dura-removal `inOriginal` input
should be the masked N4 image or an unmasked, grid-matched N4 image, because
the historic SPECTRE helper returns both a stripped image and a transferred
original image.

`run_Toads_jist.m` explicitly passes atlas scale `2.0`, smoothing `0.2`,
maximum change `1e-4`, 50 iterations, rigid registration initialization, and
`(18,6)` connectivity. Other relevant flags are `inCorrect=true`,
`inOutput2=false`, and `inOutput3=false`. These values are fixed in the
selected wrapper and must be checked against `MedicAlgorithmToads08` before
they become Python defaults. `run_RemoveDura_jist.m` passes curvature `0.5`,
level `1.0`, alignment `0.5`, pressure `0.5`, target `0.3`, mask `2.0`,
secondary level `0.0`, maximum `5`, WM threshold `0.05`, and GM threshold
`0.02`.

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
