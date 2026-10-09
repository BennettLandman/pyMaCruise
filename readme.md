# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-09

The corrected seven-stage brain-resolution synthetic base MaCRUISE run completed and resumes cleanly. Both TOADS passes naturally converged (11 and 43 iterations); an independent audit verified all 46 output hashes, all 31 finite RAS NIfTIs, and RAS-millimetre mesh headers. Peak sampled RSS was 8.50 GB for the run and 0.121 GB for resume. See the [corrected base workflow report](docs/reports/base-workflow-retained-center-2026-10-09.md).

The retained-registration-center fix has exact held-fixed Java comparisons for both TOADS passes (0/6,092,928 differing voxels; minimum label Dice 1.0). The integrated second-pass input also matched Java exactly at its natural 43-iteration stop. The current-run Dura replay matched both topology maps exactly; its brain mask differed at 3 voxels (Dice 0.9999986), within the agreed segmentation and discrete boundary-position gates. The current integrated TGDM core inputs also pass Java stage comparisons: tissue-sign is exact, and all three delivered surfaces meet Dice, mean-distance, and Hausdorff limits. See the [integrated TOADS replay](docs/reports/toads-integrated-second-source-replay-2026-10-09.md), [Dura comparison](docs/reports/dura-integrated-workflow-java-comparison-2026-10-09.md), and [TGDM comparison](docs/reports/tgdm-integrated-base-source-comparison-2026-10-09.md). The full suite on `b785bff` passed **1,503 tests, 0 failed, 0 skipped, and 11 expected RAS warnings** in 1,764.57 s, with 13.787 GB peak RSS. The coarse selected TOADS source-method estimate remains **3,511–4,911 of 13,651 lines (26–36%)**; this milestone adds integration evidence, not newly mapped Java method lines. See the [full-suite report](docs/reports/full-pytest-2026-10-09-retained-center.md) and [source-fix report](docs/reports/toads-retained-registration-center-source-fix-2026-10-09.md).

The corrected-code folded synthetic plus run also completed and resumed all nine stages, naturally converging at LesionTOADS iteration 6. Held-fixed checks on that exact run show selected Java tissue fusion, topology/refinement, and ACE are exact; GVF, TGDM, and all three JIST/Python meshes meet their source comparison criteria. The run stayed below 0.8 GB RSS, though the registered LesionTOADS loop has limited progress visibility and took about 51 minutes. These are selected-stage checks, not whole-workflow equivalence or real-MRI validation. Remaining work includes broader TOADS/LesionTOADS source coverage, plus composite parity, standalone CLI coverage, paper annotations, and real-MRI validation; the supplied human T1 still needs its matching brain mask. See the [corrected plus report](docs/reports/plus-workflow-retained-center-2026-10-09.md) and the [port plan](plan.md).

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
