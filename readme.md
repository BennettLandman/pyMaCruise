# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-11

The restartable human BASE workflow completes on hash-verified T1, brain mask, and MaSeg inputs under 16 GB; TOADS converges naturally and selected transfer/crop/surface stages match their Java references. Overall human quality remains unvalidated: tissue support is 82.4% GM and 88.4% WM (not Dice), each surface has 136 components, and the supplied Singularity output directory is empty, so the requested whole-run Dice and surface tolerances cannot yet be checked. See the [human BASE report](docs/reports/in-vivo-base-rasgrid-run-attempt-02-2026-10-11.md) and [TGDM island lineage](docs/reports/human-tgdm-island-lineage-2026-10-11.md).

On the synthetic plus case, selected Java comparisons now cover iteration-5→6 thinning/statistics/membership, iteration-6 growth and relations, registered shape sampling, and the full iteration-7 thinning state. The focused handoff passed byte-for-byte in **200.71 s** under a **4 GB** watchdog (0.936 GB peak); a discovered uppercase descriptor-name bug in lesion scoring is fixed and its regression passed. A corrected-type rerun changed 1,306 segmentation labels and 3,818 relation-map voxels but left the ratio at **0.131597** versus **0.0001**. The updated ribbon audit still fails (0% GM support, 59.77% WM support, two components per mesh, median thickness 0.373 mm). The prior objective remains Python-derived, and full live Java-loop equivalence remains open. See the [corrected P41z rerun](docs/reports/p41z-case-normalization-rerun-2026-10-11.md). Latest full suite: **1,531 passed** on `9cade54` (13.450 GB peak); no full suite has run on newer changes. Selected TOADS source-method mapping remains **3,511–4,911 of 13,651 lines (26–36%)**. See the [full-suite report](docs/reports/full-pytest-2026-10-11-toads-skeleton-carry.md), [P41z core replay](docs/reports/p41z-iteration5-to6-core-replay-2026-10-11.md), [P41z iteration handoff](docs/reports/p41z-iteration6-membership-carry-2026-10-11.md), and [detailed plan](plan.md).

## Scope and target

- **Base workflow:** T1-weighted MRI, MaSeg labels, and a user-provided brain mask on the T1 grid.
- **Plus workflow:** T1, MaSeg, matching brain mask, and FLAIR, following the selected `run_maCRUISEplus_v1` profile in [docs/plus-baseline.md](docs/plus-baseline.md).
- **Reference source:** MATLAB MaCRUISE V3 and the supplied VandyCRUISE algorithms. Numbered historical runners are bug-fix versions; versions before V3 contain known bugs.
- **Goal:** run the selected workflows without MATLAB or a MIPAV/JIST JVM runtime. ANTsPy/ANTs is the selected external backend for ANTs operations.
- **Excluded:** visualization porting and FreeSurfer requirements. A post-port discussion will select modern Python visualization tools.

## Implementation principles

1. Use existing Python libraries for general image I/O, array arithmetic, registration, filtering, morphology, and mesh operations when their behavior fits the source; port source-specific behavior with direct evidence where it does not.
2. Preserve useful MaCRUISE outputs, labels, data types, transforms, and source parameters. Treat image calculation as a generic geometry-aware N-dimensional array operator.
3. Save image inputs, fixtures, intermediates, and outputs as compressed NIfTI (`.nii.gz`). Report spatial outputs in RAS; warn and record reorientation of non-RAS inputs.
4. Build and validate on deterministic synthetic data before using the supplied human scan. Keep synthetic fixtures and generators in `tests/data/synthetic/`.
5. Report source-equivalence evidence, Python integration regressions, and real-data validation as distinct evidence categories. Do not infer whole-workflow equivalence from component matches.

The runtime stack is NiBabel, NumPy, and SciPy; pytest, build, and scikit-image are development dependencies, and `antspyx` is optional for ANTs-backed operations. The license/attribution inventory and follow-up TTO work are in [license_info.md](license_info.md). All atlas images needed at runtime are bundled in this repository as regular Git files.

## Inputs and outputs

The command-line input, assumption, and output contract is documented in [docs/cli-help.md](docs/cli-help.md). Coordinate and NIfTI rules are in [docs/image-coordinate-conventions.md](docs/image-coordinate-conventions.md). Numerical acceptance thresholds are in [docs/validation-tolerances.md](docs/validation-tolerances.md).
