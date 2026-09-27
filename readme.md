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
2. Preserve image geometry, voxel types, labels, coordinate transforms, and useful intermediate/output products. Do not reduce volumes to 2D display images or silently resample mismatched inputs.
3. Use standard Python libraries for generic operations. Image arithmetic such as the old Image Calculator step is a generic N-dimensional array operation with geometry checks and metadata preservation.
4. Use ANTsPy for ANTs operations. Replace legacy NiftyReg/Aladin MNI registration with ANTsPy only after transform direction, interpolation, and output alignment are validated against a reference case.
5. Keep TOADS, topology correction, ACE, GVF, and CRUISE/TGDM faithful to VandyCRUISE. Generic morphology or active-contour functions are not assumed to be equivalent substitutes.
6. Keep visualization and FreeSurfer comparison utilities optional and outside the core segmentation path.
7. Record source version, parameters, and transformations in run metadata for reproducibility.

## Proposed Python stack

- **NiBabel:** NIfTI and surface-file I/O, image headers/affines, and FreeSurfer formats.
- **NumPy:** volume arrays and generic voxel-wise operators.
- **SciPy ndimage:** standard filtering, connected components, interpolation, and morphology.
- **ANTsPy (`antspyx`):** N4 bias correction and ANTs registration operations.
- **scikit-image:** standard image processing only where its behavior matches the source step.
- **VTK or PyVista:** VTK surface I/O, mesh operations, and optional visualization.
- **pytest:** focused unit, stage, and end-to-end fixtures when implementation validation begins.

Pin compatible package versions after establishing the implementation environment.

## Plan and progress

See [plan.md](plan.md) for the authoritative step-by-step porting table, per-step validation criteria, synthetic-first development gates, and start/validation dates. Status key in that table: **Done**, **Started**, **Blocked**, **Not started**.

Current progress: `run_maCRUISE_v3` at MaCRUISE SVN revision 9360 is selected as the corrected base reference. The project owner clarified that numbered runners are bug-fix versions and versions below V3 had bugs; see [docs/source-baseline.md](docs/source-baseline.md). The full MaSeg-informed plus contract is documented in [docs/plus-baseline.md](docs/plus-baseline.md), and active base-stage source mappings are in [docs/source-map.md](docs/source-map.md). The package scaffold and generic volume operator exist; synthetic fixtures, algorithm ports, and their validation have not yet been completed.

### 2026-09-27 roadmap update

- Added [plan.md](plan.md) with 47 ordered, independently validatable porting steps, validation gates, and date fields.
- Split base MaCRUISE and MaCRUISEplus into explicit source-mapped implementation tracks.
- Added a synthetic-data readiness gate before any real MRI integration, a committed-fixture policy, stage-wise numerical equivalence criteria, and final real-data acceptance phases.
- Completed roadmap step P02: selected `run_maCRUISE_v3` at SVN revision 9360 as the canonical base behavior and recorded the variant comparison and effective parameters in `docs/source-baseline.md`.
- Clarified the version history: numbered runner versions are bug fixes; versions below V3 had bugs. V3 is the corrected baseline, and its observed topology/boundary behavior is part of that contract rather than an optional algorithm generation.
- Completed roadmap step P03: selected the full MaSeg-informed T1 + MaSeg + FLAIR plus workflow using `run_maCRUISEplus_v1`; documented its settings and the distinct deferred Origin/no-MaSeg variants in `docs/plus-baseline.md`.
- Completed roadmap step P04: traced the selected V3 pipeline and outputs to Python libraries, VandyCRUISE algorithm families, or explicit source gaps in `docs/source-map.md`. Custom JIST adapter sources needed to reproduce several mandatory stages are not present in the supplied source trees and are listed for P05 follow-up.
- Existing package scaffold and generic volume operator are marked started, not validated; no tests have been run yet.

## Validation plan

Validation is required at image, algorithm-stage, surface, and full-pipeline levels. Visual plausibility alone is insufficient.

### Input and image integrity

- Confirm NIfTI shape, affine, voxel spacing, orientation, units, header fields, and data type are read and written consistently.
- Reject or explicitly resolve mismatched image/mask geometry; never silently crop to smaller array dimensions.
- Check that binary masks remain binary and label maps retain their label set unless a stage explicitly changes it.
- Verify image arithmetic uses intended mask/value scaling and does not unintentionally clamp or quantize intensities.

### Stage-level equivalence

- Compare each ported stage with selected MATLAB/JIST reference output on the same input and parameters.
- Use exact equality for deterministic label/mask transforms where appropriate; document numerical tolerances for floating-point operations.
- Report Dice/Jaccard overlap, voxel-wise differences, volume/label counts, and spatial alignment for segmentation outputs.
- Compare registration transforms and resampled outputs, including transform direction, interpolation, and reference grid.
- Record parameters, software versions, and transform files for each comparison.

### Topology, surfaces, and thickness

- Check connected components, holes/Euler characteristic, and required foreground/background connectivity after topology-sensitive stages.
- Compare surface vertex/face counts, connected components, normals/orientation, coordinate bounds, and physical-space alignment.
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
- **FreeSurfer:** not required for core segmentation; used by optional FreeView rendering and comparisons with subject-specific FreeSurfer surfaces.
- **Reference data needed:** one de-identified T1 + MaSeg case, expected outputs/intermediates, and exact configuration. Add FLAIR and plus-workflow outputs for MaCRUISEplus.

## Progress log

### 2026-09-27

- Created the initial plan before implementation.
- Recorded decisions to use VandyCRUISE, use ANTsPy, and treat image arithmetic as a generic operator.
- Confirmed that source checkouts are alongside `pymacruise`; MATLAB source is in the separate `deleteme/masimatlab` tree.
- Added minimal package metadata and a geometry-checked generic NIfTI volume operator in `src/pymacruise/volume_ops.py`.
- The volume operator supports add, subtract, multiply, divide, minimum, and maximum; it rejects mismatched grids rather than implicitly resampling.
- No output-equivalence validation has been completed yet.
