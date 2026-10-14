# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-14

The port has runnable restartable BASE and PLUS workflows, selected-stage Java comparisons, and committed synthetic fixtures. Acceptance means reproducing the selected original Java/MATLAB workflows within the agreed tolerances, including their limitations. Anatomical quality audits are diagnostic; improving on the original is outside the port scope. See the [owner clarification and comparison policy](docs/validation-tolerances.md#port-acceptance-objective--owner-clarification-2026-10-14). Whole-workflow equivalence remains unverified:

- **Synthetic:** BASE completes and resumes with natural TOADS convergence. The 172×216×164 PLUS case converges and passes selected Java TGDM comparisons, but its anatomical diagnostics flag GM support, nesting, mesh topology, and thickness; these flags do not by themselves fail port acceptance. See [PLUS terminal evidence](docs/reports/plus-atlas-consistent-corrected-terminal-2026-10-12.md).
- **Real human BASE:** Case 100249 naturally converges and its RAS and LPS-compatible TOADS routes match selected Java exactly. A fresh case-100260 no-regrid BASE run with corrected Dura settings completes all eight stages, verifies all 68 stage outputs on resume, and naturally converges in TOADS at iterations 21/18. The corrected second TOADS pass matches selected Java exactly on the same Dura handoff (all 11 labels, zero differing voxels). The independent comparison against supplied historical outputs still fails the 0.97 per-label and surface-distance limits; the output tree has no run manifest, so its source settings are unverified. Current anatomy diagnostics also show incomplete ribbon support and fragmented meshes. This is human integration plus selected-stage source evidence, not whole-workflow equivalence. See the [fresh run and audit](docs/reports/in-vivo-base-v3-dura002-no-regrid-2026-10-14.md), [human Dura comparison](docs/reports/human-dura-matched-java-comparison-2026-10-14.md), and [fresh MATLAB refinement comparison](docs/reports/matlab-refinement-direct-source-comparison-2026-10-14.md).
- **Real human PLUS:** Case 100260 completes and resumes all nine stages with both geometry modes and natural iteration-13 LesionTOADS convergence. A fresh axis-only run reproduces the previous run's six TGDM fields exactly; corrected crop replay also matches all six saved TGDM products bit-for-bit. Tissue fusion and ACE match selected Java exactly, GVF is within bounds, and all five TGDM surfaces pass the agreed distance gates. A complete independent Java/MATLAB PLUS reference is still unavailable. See the [workflow audit](docs/reports/in-vivo-100260-plus-full-run-2026-10-14.md) and [human held-fixed comparisons](docs/reports/in-vivo-100260-plus-held-fixed-java-comparison-2026-10-14.md).
- **Validation:** After the selected V3 Dura WM/GM regularization correction, the full suite passed **1,577 tests, 4 skipped, 0 failed** on implementation commit `9ee3ee92965a7b71f8107086dc51334f9ffd6935` in 719.49 seconds; peak sampled RSS was 15.100 GB under a 16 GB watchdog. See the [suite report](docs/reports/full-pytest-2026-10-14-dura-v3-regularization.md). These checks separate synthetic execution, selected-source equivalence, and human anatomy; they do not establish whole-workflow equivalence.

Next work is to complete matched BASE and PLUS source comparisons, investigate Python-versus-original differences, and close remaining standalone-tool/documentation tasks. Use human ribbon-quality diagnostics to localize discrepancies rather than require better anatomy than the original. The [plan](plan.md) tracks steps and evidence; detailed comparisons are in [docs/reports](docs/reports/). Visualization remains excluded from the port, and no FreeSurfer runtime is required.

Primary validation uses our own controlled runs of the supplied original Java/MATLAB implementations. Supplied precomputed outputs are secondary comparisons for identifying areas of concern; their uncertain history does not replace fresh source execution. See the [reference evidence hierarchy](docs/validation-tolerances.md#reference-evidence-hierarchy--owner-direction-2026-10-14).

Fresh execution of the original MATLAB refinement routine identified a Python `tabulate` interpretation defect, now corrected. Small, brain-resolution, and refreshed case-100260 human refinements match MATLAB exactly on held-fixed inputs. The corrected-Dura full BASE run and its distinct historical-output discrepancies are detailed in the [run report](docs/reports/in-vivo-base-v3-dura002-no-regrid-2026-10-14.md).

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
