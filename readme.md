# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-06

The synthetic base workflow completes all six stages at the 172×216×164 input grid, with both TOADS passes converged (iterations 39 and 25), 30 finite RAS `.nii.gz` outputs, and verified resume. The raw-input plus CLI completes all nine stages on committed 25³ and 36×45×34 folded synthetic cases. P41v confirms the larger case reaches natural LesionTOADS convergence at iteration 5 with default tolerance `1e-4`, writes terminal classification/surfaces/thickness outputs, and resumes; P42a5 repeats this from an isolated installed wheel outside the checkout and verifies all nine stages and outputs. P41w's 64×80×64 plus attempt has not completed: it reached the registered LesionTOADS branch, sampled 25.825 GB peak RSS, and exited `-9` below the configured 32 GB watchdog limit. The cause is unconfirmed; reducing memory in that branch is now the next integration milestone. These are synthetic Python integrations, not whole-workflow source equivalence or real-MRI validation. In the compact Vandy synthetic comparison, Java and Python match 19 memberships and rigid poses bit-for-bit through natural convergence at iteration 19, with identical final labels (Dice 1.0 per label); the sparse TOADS queue also matches source-derived traces of 343 thinning and 480 growing events. The latest full suite on `09a7bd0` passed **1,219 tests, 0 failed, 0 skipped, 11 expected RAS warnings in 1,699.06s**, peak RSS **13.052 GB** under the 24 GB watchdog. Selected TOADS source behavior remains estimated at **3,511–4,911 of 13,651 lines (about 26–36%)**; this resource investigation adds no new TOADS source lines. Major open work is complete source-behavior coverage, base/plus whole-workflow source comparison and numerical acceptance, plus resource safety at larger grids, and real-MRI validation. No human MRI was used. See the [P41v natural-convergence report](docs/reports/plus-natural-convergence-2026-10-06.md), [P42a5 installed-wheel report](docs/reports/installed-wheel-plus-natural-convergence-2026-10-06.md), [P41w/P41x resource investigation](docs/reports/plus-scaled-grid-resource-2026-10-06.md), [compact TOADS Java comparison](docs/reports/toads-jama-lu-compact-vandy-iteration-trace-2026-10-06.md), [compiled base report](docs/reports/toads-compiled-queue-base-run-2026-10-06.md), and [detailed plan](plan.md).
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

The core stack is NiBabel, NumPy, SciPy, scikit-image, and pytest; `antspyx` is optional and required for ANTs-backed operations. The license/attribution inventory and follow-up TTO work are in [license_info.md](license_info.md). All atlas images needed at runtime are bundled in this repository as regular Git files.

## Inputs and outputs

The command-line input, assumption, and output contract is documented in [docs/cli-help.md](docs/cli-help.md). Coordinate and NIfTI rules are in [docs/image-coordinate-conventions.md](docs/image-coordinate-conventions.md). Numerical acceptance thresholds are in [docs/validation-tolerances.md](docs/validation-tolerances.md).
