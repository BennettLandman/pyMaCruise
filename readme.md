# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-11

The corrected human base tail completed from the hash-verified TOADS checkpoint with the cleaned MaSeg CSF boundary active in ACE mode 2. Preparation, CRUISE core, refinement, and surface/thickness stages completed; a same-settings resume reused all four stages. The human result still misses the requested anatomy target: GM/WM containment is **82.38%/88.49%**, not Dice, and each surface has 143 connected components. Direct Java surface extraction on these fixed fields matches Python counts and passes the distance limits, but does not establish full-workflow equivalence. The earlier human metrics came from a mode-2 run missing its boundary input and are marked preliminary. See the [corrected human-tail report](docs/reports/in-vivo-base-tail-with-v3-boundary-2026-10-11.md), [boundary integration evidence](docs/reports/base-ace-boundary-integration-2026-10-11.md), and [detailed plan](plan.md).

The latest full suite on tested implementation commit `f3b11c6` reports **1,512 passed, 19 failed, 0 skipped, 11 warnings** in 781.81 s, with 13.535 GB peak sampled RSS under a 16 GB watchdog. All failures are TOADS thinning/iteration source-reference checks; two reproduce individually. They follow the TOADS constructor-border-state change at `c3c0478` and are under source/harness audit; the suite does not pass. See the [full-suite report](docs/reports/full-pytest-2026-10-11-v3-boundary.md). The selected TOADS source-method estimate remains **3,511–4,911 of 13,651 (26–36%)**, with zero mapped Java lines added by boundary integration. Remaining work includes human ribbon-quality acceptance, Dura parameter interpretation, historical full-pipeline comparison, plus-workflow quality/parity, broader TOADS/LesionTOADS behavior, standalone CLIs, and paper annotations. Synthetic, fixed-component source comparison, and human integration evidence remain separate.

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
