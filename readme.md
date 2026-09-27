# PyMaCRUISE

Python-native port of the MaCRUISE cortical segmentation and surface-generation workflow. The goal is to preserve the behavior and useful outputs of the existing MATLAB/JIST/MIPAV workflow while reusing established Python imaging libraries wherever they fit. We will port MaCRUISE-specific algorithms and workflow logic that do not have a suitable reusable implementation; routine file I/O, image operations, registration, and mesh handling should use existing libraries or supported tool wrappers.

This README is the project plan and progress record. Update the checklist and decisions here as implementation and validation proceed.

## Scope and target

- **Primary workflow:** MaCRUISE with a T1-weighted MRI and existing multi-atlas segmentation (`MaSeg`) as inputs.
- **Follow-on workflow:** MaCRUISEplus with T1, MaSeg, and FLAIR inputs, including lesion-aware processing.
- **Reference MATLAB source:** `../../deleteme/masimatlab/trunk/xnatspiders/matlab/MaCRUISE/`.
- **Selected algorithm source:** `../toads-cruise-scm-2026-09-26/VandyCRUISE/src/`.
- **Registration/bias-field toolkit:** ANTsPy (`antspyx`) where it can reproduce the reference operation.
- **Final runtime goal:** Python workflow and Python-accessible algorithm implementations, without MATLAB or a MIPAV/JIST JVM runtime in the completed pipeline. Native library backends such as ITK/ANTs are acceptable.

Start with the base T1 + MaSeg workflow. Add MaCRUISEplus after the base path and its validation are stable.

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

Status key: `[x]` complete, `[~]` in progress, `[ ]` not started.

### 0. Project setup

- [x] Create this plan and progress record.
- [x] Establish the initial Python package structure and dependency metadata.

### 1. Establish contracts and source mapping

- [ ] Identify the canonical MATLAB entry point and pipeline variant for the first release.
- [ ] Map active MATLAB steps, parameters, inputs/outputs, and Java classes to VandyCRUISE algorithms or Python/library replacements.
- [ ] Document JIST/MIPAV wrapper behavior that must be re-created in Python.
- [ ] Inventory required atlas/configuration assets and external dependencies; identify what is bundled.
- [ ] Define input contracts (T1, MaSeg, geometry, labels) and versioned outputs/intermediates.

### 2. Create Python package and workflow foundation

- [ ] Add package metadata, supported Python version, dependency pins, and a command-line entry point.
- [ ] Add typed configuration for inputs, output directory, parameters, and optional stages.
- [ ] Implement geometry-aware volume I/O and validation of shape, affine, orientation, spacing, and data type.
- [ ] Implement run manifests, structured logging, stage output directories, and clear failure reporting.
- [x] Implement generic voxel operators with explicit geometry and dtype rules.

### 3. Port the base MaCRUISE stages

- [ ] Raw input preparation and orientation/metadata handling.
- [ ] Cortical label separation and mask preparation.
- [ ] TOADS preprocessing/segmentation and transfer back to input space; use ANTsPy for N4 where appropriate.
- [ ] Crop/uncrop, atlas/TOADS combination, and multi-atlas boundary/CSF preparation.
- [ ] Topology correction.
- [ ] ACE tissue enhancement and skeleton generation.
- [ ] Gradient Vector Flow (GVF).
- [ ] CRUISE/TGDM evolution of inner, central, and outer surfaces/level sets.
- [ ] Mesh extraction, transfer to full/original space, cortical thickness, and MNI outputs.
- [ ] Preserve useful intermediates and established output naming where practical.

### 4. Add the MaCRUISEplus path

- [ ] FLAIR-to-T1 rigid registration with ANTsPy and reference-compatible settings.
- [ ] Lesion TOADS and associated masking/cleanup.
- [ ] Plus-specific boundary, ACE, topology, and CRUISE/TGDM parameter differences.
- [ ] Validate plus outputs separately from the base workflow.

### 5. Optional outputs and usability

- [ ] Python surface/volume visualization and report generation to replace optional FreeView screenshots.
- [ ] Optional comparison with FreeSurfer subject surfaces and regional thickness results.
- [ ] Batch/cluster execution support after the single-subject workflow is stable.
- [ ] Installation, input, command, output, and troubleshooting documentation.

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
