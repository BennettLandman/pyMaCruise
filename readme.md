# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-14

The port has runnable restartable BASE and PLUS workflows, selected-stage Java comparisons, and committed synthetic fixtures. Current evidence is not a full acceptance result:

- **Synthetic:** BASE completes and resumes with natural TOADS convergence. The 172×216×164 PLUS case converges and passes selected Java TGDM comparisons, but misses GM support, nesting, mesh-topology, and thickness checks. See [PLUS terminal evidence](docs/reports/plus-atlas-consistent-corrected-terminal-2026-10-12.md).
- **Real human BASE:** Case 100249 naturally converges and its RAS and LPS-compatible TOADS routes match selected Java exactly. Case 100260 completes all eight stages with 68/68 hashes and 51/51 finite RAS NIfTIs; its second TOADS pass reaches the source 50-iteration cap, with exact Java classification parity at that cap. Its ribbon support and surface connectivity checks fail. See the [100260 run audit](docs/reports/in-vivo-100260-base-run-audit-2026-10-14.md) and [Java comparison](docs/reports/in-vivo-100260-second-toads-java-comparison-2026-10-14.md).
- **Validation:** The full suite passed **1,573 tests, 4 skipped, 0 failed** on implementation commit `0981b8d153b700615450ffa64844b2404d584137` in 708.85 seconds; peak sampled RSS was 14.857 GB under a 16 GB watchdog. See the [suite report](docs/reports/full-pytest-2026-10-14-toads-lps-axis-compatibility.md). These checks separate synthetic execution, selected-source equivalence, and human anatomy; they do not establish whole-workflow equivalence.

Next work is to resolve or characterize the human ribbon-quality failures, finish PLUS integration and source checks, and close remaining standalone-tool/documentation tasks. The [plan](plan.md) tracks steps and evidence; detailed comparisons are in [docs/reports](docs/reports/). Visualization remains excluded from the port, and no FreeSurfer runtime is required.

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
