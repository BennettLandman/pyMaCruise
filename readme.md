# PyMaCRUISE

Python-native port of the MaCRUISE cortical segmentation and surface-generation workflow. The goal is to preserve the behavior and useful outputs of the existing MATLAB/JIST/MIPAV workflow while reusing established Python imaging libraries wherever they fit. We will port MaCRUISE-specific algorithms and workflow logic that do not have a suitable reusable implementation; routine file I/O, image operations, registration, and mesh handling should use existing libraries or supported tool wrappers.

The detailed, step-by-step porting roadmap and acceptance table are in [plan.md](plan.md). This README is the project entry point and concise progress/decision summary; update it alongside the roadmap as implementation and validation proceed.

## Scope and target

- **Primary workflow:** MaCRUISE with a T1-weighted MRI and existing multi-atlas segmentation (`MaSeg`) as inputs.
- **Follow-on workflow:** MaCRUISEplus with T1, MaSeg, and FLAIR inputs, including lesion-aware processing.
- **Reference MATLAB source:** `../../deleteme/masimatlab/trunk/xnatspiders/matlab/MaCRUISE/`.
- **Selected algorithm source:** `../toads-cruise-scm-2026-09-26/VandyCRUISE/src/`.
- **Registration/bias-field toolkit:** ANTsPy (`antspyx`) where it can reproduce the reference operation.
- **Final runtime goal:** Python workflow and Python-accessible algorithm implementations, without MATLAB or a MIPAV/JIST JVM runtime in the completed pipeline. Native library backends such as ITK/ANTs are acceptable.

Start with the base T1 + MaSeg workflow. The selected follow-on plus contract is full MaSeg-informed T1 + MaSeg + FLAIR using the `_v1` wrapper profile; details and deferred alternatives are in [docs/plus-baseline.md](docs/plus-baseline.md). Add its implementation after the base path and validation are stable.

## Porting principles

1. Treat the MATLAB files as the workflow and parameter reference; treat VandyCRUISE as the selected algorithm source.
2. Preserve image geometry, voxel types, labels, coordinate transforms, and useful intermediate/output products. Report all outputs in RAS; reorient non-RAS inputs to RAS with a visible warning and recorded transform. Do not reduce volumes to 2D display images or silently resample mismatched inputs.
3. Use standard Python libraries for generic operations. Image arithmetic such as the old Image Calculator step is a generic N-dimensional array operation with geometry checks and metadata preservation.
4. Use ANTsPy for ANTs operations. Replace legacy NiftyReg/Aladin MNI registration with ANTsPy only after transform direction, interpolation, and output alignment are validated against a reference case.
5. Keep TOADS, topology correction, ACE, GVF, and CRUISE/TGDM faithful to VandyCRUISE. Generic morphology or active-contour functions are not assumed to be equivalent substitutes.
6. Record source version, parameters, and transformations in run metadata for reproducibility. Visualization is outside this port; discuss modern Python tools after acceptance at P47.

## Proposed Python stack

- **NiBabel:** NIfTI I/O, image headers/affines, and geometry metadata.
- **NumPy:** volume arrays and generic voxel-wise operators.
- **SciPy ndimage:** standard filtering, connected components, interpolation, and morphology.
- **ANTsPy (`antspyx`):** N4 bias correction and ANTs registration operations.
- **scikit-image:** standard image processing only where its behavior matches the source step.
- **Scientific surfaces:** select a mesh geometry/format library only if required by P30–P31; no renderer or viewer dependency is in scope.
- **pytest:** focused unit, stage, and end-to-end fixtures.

Core dependency compatibility bounds and the optional ANTs extra are declared in `pyproject.toml`; the current validated environment versions are recorded in [docs/assets-and-dependencies.md](docs/assets-and-dependencies.md).

PyMaCRUISE has no FreeSurfer dependency or runtime requirement: users will not need FreeSurfer installed, licensed, configured, or provide FreeSurfer data or outputs. Mentions in the source notes describe historical MaCRUISE behavior only.

All image data stored by PyMaCRUISE—including bundled image assets, synthetic fixtures, cached intermediates, and outputs—uses compressed NIfTI (`.nii.gz`). Bundled atlas image payloads have been normalized to RAS; the source and conversion audit is under `data/atlases/`. Non-image metadata/configuration files keep their appropriate formats. Some historical source geometry metadata remains uncertain and is clearly flagged in the audit and P05c.

## Plan and progress

See [plan.md](plan.md) for the authoritative step-by-step porting table, per-step validation criteria, synthetic-first development gates, and start/validation dates. Status key in that table: **Done**, **Started**, **Blocked**, **Not started**.

Current progress: `run_maCRUISE_v3` at MaCRUISE SVN revision 9360 is selected as the corrected base reference. The project owner clarified that numbered runners are bug-fix versions and versions below V3 had bugs; see [docs/source-baseline.md](docs/source-baseline.md). The full MaSeg-informed plus contract is documented in [docs/plus-baseline.md](docs/plus-baseline.md), active base-stage mappings are in [docs/source-map.md](docs/source-map.md), source/assets/runtime requirements are in [docs/assets-and-dependencies.md](docs/assets-and-dependencies.md), working image/coordinate conventions are in [docs/image-coordinate-conventions.md](docs/image-coordinate-conventions.md), and the input/output contract is in [docs/io-contract.md](docs/io-contract.md). Supplied TOADS/CRUISE, SPECTRE, and lesion-TOADS runtime assets are bundled under `data/atlases/` as regular Git files; all 35 atlas NIfTIs have exact voxel-value matches and RAS voxel axes, with source geometry caveats recorded. The license/attribution inventory and TTO handoff notes are in [license_info.md](license_info.md); preserved upstream notice files are under `licenses/upstream/`. The MNI reference/mask are absent from the supplied source trees. Synthetic fixtures and generic infrastructure are in place; MaCRUISE algorithm stages and equivalence checks remain to be ported.

The currently implemented CLI entry point supports the implemented stages:

```text
pymacruise prepare-raw --t1 T1.nii.gz --maseg MaSeg.nii.gz --output-dir Raw
pymacruise cortical-separation --maseg Seg.nii.gz --output-dir CorticalSeparation
pymacruise calculate --input-a A.nii.gz --input-b B.nii.gz --operation ADD --output Sum.nii.gz
```

### 2026-09-27 roadmap update

- Added [plan.md](plan.md) with ordered, independently validatable porting steps, validation gates, and date fields.
- Split base MaCRUISE and MaCRUISEplus into explicit source-mapped implementation tracks.
- Added a synthetic-data readiness gate before any real MRI integration, a committed-fixture policy, stage-wise numerical equivalence criteria, and final real-data acceptance phases.
- Completed roadmap step P02: selected `run_maCRUISE_v3` at SVN revision 9360 as the canonical base behavior and recorded the variant comparison and effective parameters in `docs/source-baseline.md`.
- Clarified the version history: numbered runner versions are bug fixes; versions below V3 had bugs. V3 is the corrected baseline, and its observed topology/boundary behavior is part of that contract rather than an optional algorithm generation.
- Completed roadmap step P03: selected the full MaSeg-informed T1 + MaSeg + FLAIR plus workflow using `run_maCRUISEplus_v1`; documented its settings and the distinct deferred Origin/no-MaSeg variants in `docs/plus-baseline.md`.
- Completed roadmap step P04: traced the selected V3 pipeline and outputs to Python libraries, VandyCRUISE algorithm families, and the MaCRUISE-specific adapters in `toads-cruise/src` in `docs/source-map.md`.
- Completed roadmap step P05: cataloged the source trees, CVS/SVN version evidence, atlas payloads, available notices, absent MNI references, external-tool availability, and dependency gaps in `docs/assets-and-dependencies.md`.
- Completed roadmap step P05a: copied all supplied atlas trees required by the selected base/plus wrapper paths into `data/atlases/`, excluded CVS working-copy metadata, stored every payload as a regular Git file, and recorded SHA-256 entries. The MNI reference/mask remain unavailable in the supplied trees.
- Completed roadmap step P05b: preserved all 70 files from the supplied MIPAV/JIST license directories under `licenses/upstream/` and documented source attributions, module-comment practice, owner direction, and the later TTO licensing handoff in `license_info.md`.
- Started roadmap step P06: recorded source-backed NIfTI array, orientation-transfer, coordinate-frame, interpolation, crop, and header rules in `docs/image-coordinate-conventions.md`. Affine/axis, landmark, crop, and transform checks remain pending; P06 is not validated.
- Completed P07: added profile-aware synthetic T1/MaSeg/FLAIR validation for dimensionality, finite intensity values, integer labels, grid agreement, and RAS warning/reorientation behavior; documented the full input/output policy.
- Completed P08: the user confirmed minimum per-region Dice ≥0.97, symmetric surface Hausdorff ≤1 mm, and symmetric mean nearest-surface distance ≤0.01 mm; output-class comparison methods and synthetic metric evidence are recorded in `docs/validation-tolerances.md`.
- Completed roadmap steps P09–P10: added a seeded phantom generator and committed a compact, checksummed corpus of compressed NIfTI image fixtures plus a RAS surface mesh under `tests/data/synthetic/`.
- Visualization tools, FreeView scripting/rendering, and FreeSurfer comparisons are excluded from this port. P47 reserves a post-acceptance discussion of modern Python visualization options.
- Completed P11: added exact voxel, Dice/Jaccard, voxel error, point-cloud/triangle-surface, transform, and mesh-topology metrics plus deterministic JSON reports. The current offline synthetic suite passes 59 tests.
- Started P05c: normalized 35 user-provided atlas NIfTIs to RAS, moved them to their source atlas folders, updated descriptors, removed legacy image payloads from the working tree, and regenerated checksums. Exact voxel data checks pass for all assets; original world-coordinate certainty remains unresolved for Analyze examples and source files with incomplete MIPAV metadata.
- Started P12: bounded core/optional dependencies in `pyproject.toml`; editable installation succeeded in the project virtual environment without installing ANTsPy. A fresh-environment installation check remains.
- Completed P12–P14: bounded core/optional dependencies, added the `pymacruise` CLI entry point, implemented geometry-aware NIfTI I/O, validated typed workflow profiles, and added run manifests with dependency versions, hashes, and configuration-sensitive signatures.
- Completed P15–P17: verified generic voxel arithmetic; ported RAS-normalized `prepare_raw`; ported exact source-label MaSeg cortical separation and all four masks.
- Started P22: source-traced the default 3-D boundary algorithm and selected CSF merge rules; its thin-skeleton product is separately gated with P27.
- Completed P22: ported the source's default 3-D six-neighbor erosion boundary, overlooked-label removal, low-label contact suppression, and selected CSF merge modes 1 and 3; exact synthetic arrays and RAS image geometry pass.
- Completed P23: ported both exact V3 boundary-cleanup thresholds (`>0.01` and `>25.5`), with equality/above-threshold tests and aligned RAS image handling.
- Completed P24: ported the selected V3 Java atlas/TOADS tissue combination, including integer six-neighbor WM erosion, confidence rules, branch priority, output labels/dtypes, and brain-mask behavior.
- Completed P21: ported TOADS-to-reference transfer with nearest-neighbor interpolation, using physical affine mapping to the reference grid and RAS outputs; synthetic translations, rotations, labels, and continuous values pass.
- Completed P20: ported GM-derived crop and inverse crop in `src/pymacruise/cropping.py`. Synthetic asymmetric ROI/landmark cases validate the strict threshold, inclusive padded bounds, RAS affine-origin update, outside-image fill, grid rejection, and non-RAS warnings. Bytecode inspection found a historical JIST inverse-crop off-by-one (exclusive maximum); default behavior preserves it, and a complete inverse is available explicitly. The full suite passes 63 tests.
- Completed P25: implemented the selected VandyCRUISE downward `(6,18)` topology propagation and V3 multi-atlas WM refinement in `src/pymacruise/topology.py`. Synthetic disconnected-object, torus-handle, cavity, local-connectivity, repeatability, atlas-refinement, geometry, and warning tests pass; numerical parity with the Java queue remains for reference validation. Details are in [docs/topology-correction.md](docs/topology-correction.md).
- Started P26: tracing the MaCRUISE ACE adapter and VandyCRUISE enhancement/diffusion/skeleton code. The selected profile uses threshold 126.9, Gaussian alpha 1.8, and applies updated CSF after diffusion.
- Started P18: source-reviewed canonical TOADS preprocessing: N4, FSL BET (`-f 0.1 -B`) followed by SPECTRE, two TOADS passes with dura removal between, and exact stage outputs. BET/SPECTRE execution remains outstanding; FSL BET is unavailable in the current environment.
- Completed P19: implemented the selected N4 stage using ANTsPy with settings equivalent to the legacy ANTs command defaults; synthetic bias-field, grid, repeatability, and input-contract checks pass. ANTsPy 0.6.3 runs with restored SciPy 1.18.1 here, but its published SciPy `<1.16` constraint conflicts with SciPy wheels that fail to load on this host; clean optional-extra installation remains an environment packaging issue.
- The native boundary and topology-correction algorithm stages are now implemented; the core TOADS/SPECTRE and CRUISE/TGDM algorithms remain to be ported and validated.

### 2026-09-27 implementation update

- Completed P07 and P11–P17 with profile/input contracts, a synthetic comparison harness, package/CLI scaffolding, RAS NIfTI I/O, configuration/manifests, the generic operator, raw input preparation, and exact MaSeg cortical-separation masks.
- Updated `data/atlases/`: all 35 supplied NIfTI conversions are in the correct atlas folders, use RAS axes, match the original voxel values exactly after axis reorientation, and have updated compressed-NIfTI descriptors and checksums. P05c remains Started while unresolved legacy physical-origin metadata is reviewed.
- The user confirmed the segmentation and surface-position metric interpretations; P08 is complete.
- The synthetic suite passes 59 tests; the editable install and `pymacruise --version` / `--help` succeed. The N4 test runs with the optional ANTsPy package installed.

## Validation plan

Validation is required at image, algorithm-stage, scientific surface, and full-pipeline levels, using synthetic fixtures and measurable equivalence criteria.

### Input and image integrity

- Confirm NIfTI shape, affine, voxel spacing, RAS orientation, units, header fields, and data type are read and written consistently; non-RAS inputs must warn and normalize correctly.
- Reject or explicitly resolve mismatched image/mask geometry; never silently crop to smaller array dimensions.
- Check that binary masks remain binary and label maps retain their label set unless a stage explicitly changes it.
- Verify image arithmetic uses intended mask/value scaling and does not unintentionally clamp or quantize intensities.

### Stage-level equivalence

- Compare each ported stage with selected MATLAB/JIST reference output on the same input and parameters.
- Use exact equality for deterministic label/mask transforms where appropriate; document numerical tolerances for floating-point operations.
- Report Dice/Jaccard overlap, voxel-wise differences, volume/label counts, and spatial alignment for segmentation outputs.
- Apply the project thresholds in [docs/validation-tolerances.md](docs/validation-tolerances.md): per-region Dice ≥0.97, symmetric surface Hausdorff ≤1 mm, and mean absolute surface distance ≤0.01 mm, all surface distances in RAS+ millimetres.
- Compare registration transforms and resampled outputs, including transform direction, interpolation, and reference grid.
- Record parameters, software versions, and transform files for each comparison.

### Topology, surfaces, and thickness

- Check connected components, holes/Euler characteristic, and required foreground/background connectivity after topology-sensitive stages.
- Compare surface vertex/face counts, connected components, normals/orientation, RAS+ coordinate bounds, and physical-space alignment.
- Compare inner/central/outer surface ordering and check that surfaces do not cross where the method requires nesting.
- Compare pointwise and regional thickness distributions, summary statistics, and invalid/outlier counts.

### End-to-end acceptance

- Run base MaCRUISE and MaCRUISEplus independently on de-identified reference cases with expected outputs.
- Check all required products are present, readable, aligned, and named/documented.
- Verify rerun and partial-stage restart behavior, and actionable errors for missing inputs/dependencies.
- Keep reference data and expected outputs out of version control unless explicitly approved and appropriate to share.

## Known dependencies, gaps, and decisions

- **VandyCRUISE is the selected algorithm source.** MATLAB wrappers call JIST entry-point classes. Some wrapper/orchestration classes are not under the VandyCRUISE source path and must be reconstructed from the MATLAB workflow and VandyCRUISE algorithms; do not silently switch to the separate `toads-cruise` source tree.
- **ANTsPy selected.** Confirm exact N4 and registration settings from wrappers. MNI registration currently calls NiftyReg/Aladin and needs explicit equivalence work if moved to ANTsPy.
- **Generic image arithmetic:** implement as a generic volume operator, not by porting the supplied 2D `BufferedImage` snippet. Preserve the medical-volume grid and metadata.
- **Multi-atlas fusion:** top-level MATLAB setup references separate `masi-fusion/src` code. Base MaCRUISE accepts MaSeg as input, but fusion source may be needed if producing MaSeg is in scope.
- **MNI reference:** source expects `average305_t1_tal_lin.nii.gz` and related atlas/configuration files. Locate and version these assets before validating MNI outputs.
- **Visualization:** excluded from the port. P47 is a post-acceptance discussion of modern Python tools, not a FreeView/MATLAB visualization implementation step.
- **Reference data needed:** one de-identified T1 + MaSeg case, expected outputs/intermediates, and exact configuration. Add FLAIR and plus-workflow outputs for MaCRUISEplus.

## Progress log

### 2026-09-27

- Created the initial plan before implementation.
- Recorded decisions to use VandyCRUISE, use ANTsPy, and treat image arithmetic as a generic operator.
- Confirmed that source checkouts are alongside `pymacruise`; MATLAB source is in the separate `deleteme/masimatlab` tree.
- Added minimal package metadata and a geometry-checked generic NIfTI volume operator in `src/pymacruise/volume_ops.py`.
- The volume operator supports add, subtract, multiply, divide, minimum, and maximum; it rejects mismatched grids rather than implicitly resampling.
- No output-equivalence validation has been completed yet.
