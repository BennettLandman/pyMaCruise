# PyMaCRUISE against complete external MaCRUISE

Updated October 27, 2026. This supersedes the Python-only readiness recommendation.

For BASE commands, use the [full BASE walkthrough](full-base-walkthrough.md). PLUS now has a concrete [full walkthrough](full-plus-walkthrough.md) with bundled atlas preparation, original/Python scientific commands, stage tools and explicit assumptions. The BASE container is not a verified PLUS runtime. AtlasEM/Dura and final refinement are implemented; independent source execution and whole-workflow agreement remain unverified.

## Required outcome

The tester must execute the complete external MaCRUISE workflow and the complete corresponding PyMaCRUISE workflow from the same agreed starting inputs, then compare independently generated scientific outputs. The primary deliverable is a quantitative external-versus-Python comparison report, accompanied by both full run records and visual review.

Python execution, restart checks, selected Java harnesses, and supplied historical results are supporting checks. None replaces the fresh external end-to-end run. Do not seed external scientific stages with Python-generated intermediate arrays, or seed Python scientific stages with external intermediate arrays, and call that an independent complete-workflow comparison.

## Readiness as of this revision

The current package has a synthetic-tested BASE comparator and a synthetic-tested PLUS product comparator, but **neither fresh external end-to-end workflow has been verified here**. PLUS has synthetic and human execution evidence, bundled selected-atlas conversion, and a prepared original driver. The complete original runtime and Aladin optimizer agreement remain unavailable or unverified.

Read-only repository access is sufficient for obtaining and running published code locally. It does not supply an external container/runtime. The reference-run drivers and this protocol are included in the repository. There is no requirement for the tester to commit or push changes.

## Selected external reference

The owner selected the MaCRUISE v3.2.0 container. The walkthrough pins `masidocker/harmonizeapp:macruise` to `sha256:96c6c7a033edb86f2be7b3a60c3e8cc658691de3f670b35a4147feaf287a28a6`. Its launcher scripts were inspected to call the same complete embedded application directly, preserving its exit status; outer BIDS discovery/copy/cleanup is bypassed. Record the exact external artifact, digest, entry point, and effective workflow options. Do not infer that different V3 releases or container wrappers have identical behavior.

Two public sources identify different external distributions:

- [MASILab/MaCRUISE](https://github.com/MASILab/MaCRUISE) documents `masidocker/spiders:MaCRUISE_v3_1_0`, taking `T1.nii.gz` and `orig_target_seg.nii.gz` under `/INPUTS`, with products under `/OUTPUTS`. It documents a complete container entry point and states that source is included in the container.
- [MASI harmonizeapp](https://hub.docker.com/r/masidocker/harmonizeapp) identifies `masidocker/harmonizeapp:macruise` as MaCRUISE v3.2.0 and documents a BIDS application interface. Its examples use another tag spelling, `MaCRUISEv3.2.0`; the actual selected `macruise` tag/digest and launch scripts have now been inspected rather than inferred from those examples.

The supplied local v3.2.0 log ends with missing `SEG_refine.nii.gz`; it is not a successful reference run and belongs to a different scan from the downloaded root outputs. Historical output trees are not substitutes for a newly successful run.

The project also has selected-mask MATLAB BASE/PLUS drivers in [`source-reference-run/matlab/`](../../source-reference-run/matlab/). They have static checks but have not been executed in the full MATLAB/MIPAV/JIST runtime. They are new orchestration around original algorithms and bypass original skull stripping and presentation outputs. They must not be presented as an unchanged invocation of the external product. Clone the revision containing this protocol to obtain these drivers, and record its full Git commit in the run record.

## Define the complete-workflow boundary explicitly

The external top-level application and the Python port currently have different preprocessing contracts: external variants generate skull-stripping masks, whereas Python requires a supplied mask. The Python scope also excludes legacy visualization/FreeSurfer exports. These differences must be visible in the test report.

For an **unchanged external-product test**, launch the external entry point normally and retain its internally generated preprocessing and scientific outputs. Launch Python with its declared input contract. Compare results and report preprocessing differences as actual product differences. If masks, N4 settings, registration, or geometry differ, a final discrepancy still belongs in the report; it cannot establish a like-for-like algorithm implementation defect by itself.

For **strict selected-workflow equivalence**, the starting mask, working geometry, atlas files, preprocessing choices, and algorithm settings must match by an explicitly approved shared contract. If the external product cannot accept this contract through its supported interface, an adapted driver must be separately validated and labeled as such. Feeding a mask generated by the external run into Python can help diagnose downstream differences, but changes the independence and starting boundary; report it as a separate matched-input experiment.

Neither boundary may silently replace a missing source stage, skip a failed step, use an old cached reference, or discard a difference by relaxing thresholds after seeing results.

## Tester procedure once the reference package is ready

1. **Obtain the pinned testing kit.** Clone the published Python repository revision and obtain the specified external image/runtime. Use a writable local testing directory. Record Git commit, image digest, platform, runtime/dependency versions, memory allocation, and all launch commands. The actual Python run commands and checks are in the supplemental [Python BASE checks](python-base-checks.md); its small synthetic test is installation preflight only.
2. **Verify the inputs.** Start BASE with the three files in `source-reference-run/inputs/base-root/`; use its T1 and MaSeg as external inputs, retaining the supplied mask for Python and the declared comparison boundary. Verify `source-reference-run/SHA256SUMS`. For PLUS, use the distinct `inputs/100260/` set including FLAIR. Do not mix the two scans. Create identical input copies as required by the external interface; record their hashes.
3. **Run the external pipeline first, from scratch.** Use a fresh output directory and capture full console/runtime logs. Retain every scientific intermediate, the final refined segmentation, native inner/central/outer level sets and meshes, thickness outputs, transforms, and effective settings. Confirm the required scientific outputs actually exist, are readable, and are finite. A wrapper saying "completed" is insufficient. Record a failed external run as a failure; do not continue as though a valid oracle exists.
4. **Run Python independently, from scratch.** Use a separate fresh output directory and the explicit full settings for the selected profile. Version 0.1.1 defaults use the full selected BASE profile; two-iteration limits require explicit synthetic-test overrides. Preserve its manifest, logs, effective working inputs, scientific intermediates, and outputs. Verify hashes and image integrity, then run its same-arguments resume check. For an unchanged external-product comparison, do not replace either pipeline's internal N4, segmentation, Dura, or surface products with products from the other.
5. **Verify correspondence before computing distances.** Match input identity, physical grids, orientation, interpolation, crop/uncrop mappings, mesh coordinate frames and units, atlases, stop conditions, and output definitions. Explicitly convert known source voxel-coordinate meshes to physical RAS millimetres. Do not register final outputs to each other to improve agreement. Do not compare arrays by index just because their shapes match. Any intended axis-only geometry comparison needs a documented mapping back to the same physical frame.
6. **Produce the comparison report.** Compare each required final label and each inner/central/outer surface separately. Include thickness comparisons and upstream intermediate comparisons to locate the earliest divergence. Keep all failures in the report. Review overlays in three planes and annotate discrepant locations. Numeric comparisons are primary; screenshots support interpretation.
7. **Return both run trees and the report.** Include all input/output hashes, complete logs, manifests, commands, image/runtime identity, per-label and per-surface metrics, scalar/thickness summaries, and screenshots. Return failures and incomplete runs as well. Use the project owner's agreed transfer channel; no GitHub write permission is needed. PLUS-specific product and precondition details are in the [PLUS comparison handoff](full-plus-walkthrough.md).

## Comparison criteria

These are the repository's current source-equivalence criteria; they do not mean the current end-to-end port has passed.

| Product | Required comparison |
| --- | --- |
| Final refined segmentation | Dice at least 0.97 for every required label/region; report each label and its voxel counts. An average cannot conceal a failed label. |
| Inner, central, outer surfaces | Each symmetric mean nearest-triangle-surface distance at most 0.01 mm, and symmetric sampled Hausdorff distance at most 1 mm, in the same physical RAS frame using the repository comparator definition. |
| Deterministic intermediate labels/masks | Exact equality where the selected source contract is exact; report differing voxel counts and Dice. |
| Continuous fields and thickness | MAE, RMSE, maximum absolute error, valid/invalid counts, units and distributions. Stage-specific limits must be fixed before a pass claim; no universal thickness/scalar tolerance is currently established. |
| Transforms | Frame/direction/grid checks plus rotation and translation error; agree on stage-specific limits before acceptance. |
| Anatomy/topology | Record ribbon support, nesting, components, Euler characteristic and thickness diagnostics for both pipelines. These do not substitute for source agreement. |

A report needs separate conclusions for external-run completion, Python-run completion, comparability of the workflows, segmentation agreement, surface agreement, thickness agreement, and overall acceptance. Overall acceptance remains unverified if required outputs, mappings, provenance, thresholds, or stages are missing. A Python-only run cannot pass this protocol.

## Work required before dispatching the tester

- Execute the pinned v3.2.0 reference using the inspected complete-application launch command on a suitable machine; verify its output tree and nested Docker/Singularity runtime requirements.
- Record the full published Git commit used for testing. The selected-mask reference drivers, input packet, run-record template, output auditors, and instructions are in the repository. The BASE and PLUS comparators have synthetic checks; full external launch/output-layout execution remains unverified.
- Verify `scripts/compare_external_base_run.py` against the fresh pinned-container outputs. It follows active Python attempt paths, checks input identity/output hashes, and compares final labels, level sets, surfaces and thickness. Nine synthetic checks cover its mapping and rejection behavior. The native source voxel/frame contract, fresh full-runtime execution, detailed intermediate mapping and vertexwise thickness correspondence still need verification; the historical scripts alone do not supply those checks.
- Complete one fresh external run and one independent Python run through the proposed instructions, then exercise the report generator. Record pass/fail honestly; the test can expose disagreement and still be a successful validation experiment.
- For PLUS, use the [PLUS comparison handoff](full-plus-walkthrough.md). Prepare the owner-selected `ibsr11objV3` native/Python atlas resources with the kit, obtain a working original runtime, then validate registration, AtlasEM and LesionTOADS boundaries. Current PLUS synthetic execution alone cannot be represented as a complete reproduction of that external workflow.

The host used to prepare this document has neither Docker/Apptainer/Singularity on PATH nor a demonstrated working complete original runtime. External execution and the proposed end-to-end acceptance test have not been verified here. This is the revised test contract and readiness assessment, not a claim that a turnkey external test has been completed.

## Exploratory PLUS testing policy (2026-10-27)

The owner authorizes exploratory testing before composite source equivalence is established. Follow the current PLUS walkthrough: preserve physical FLAIR geometry and record backend/atlas/transfer assumptions. The comparator defaults to strict atlas-byte identity. Explicit `--atlas-policy report-differences` computes final-product metrics while retaining atlas mismatches in the report; it does not waive input identity, physical-grid requirements or numeric gates, and cannot close whole-port acceptance. Source mesh coordinates are explicitly declared as voxel indices or RAS millimetres.
