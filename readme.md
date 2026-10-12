# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-12

The restartable BASE workflow completes on the supplied human T1 under the 16 GB memory ceiling. Orthogonal RAS regridding is optional, recommended, and enabled by default in the base and PLUS workflow commands; `--no-regrid-inputs` preserves an oblique orthogonal RAS grid without resampling. A complete no-regrid human BASE run passed all stage checks with 7.945 GB peak RSS and natural convergence in both TOADS passes. Its level-set Dice against the downloaded run meets 0.97, while surface-position and refined-label gates still fail; the fixed-field JIST mesh comparison passes. Sheared grids remain unsupported. See the [no-regrid human report](docs/reports/in-vivo-base-oblique-no-regrid-2026-10-12.md). The latest full suite passed **1,561 tests, 4 skipped, 0 failed**, with 11 expected RAS warnings, in 829.51 seconds on commit `452cc573`; peak sampled RSS was 12.018 GB under a 16 GB watchdog. See the [full-suite report](docs/reports/full-pytest-2026-10-12-oblique-optout.md) and [regridding behavior/evidence](docs/reports/orthogonal-ras-regridding-focused-2026-10-11.md).

Real-human BASE processing and selected Java stage comparisons are available, but whole-workflow source agreement and ribbon acceptance remain unresolved. Native-grid, voxel-index level-set Dice is 0.971–0.972 versus 0.943–0.957 in the previous regridded comparison; the source level-set headers do not preserve the raw oblique affine, so this is not a shared-affine comparison. Refined-label Dice and surface-position limits still fail. A selected JIST replay confirms the source N4 RAS-to-LPS x/y axis conversion (correlation 0.999998), but the physical transfer of saved TOADS/Dura products remains unresolved; intermediate stage metrics stay diagnostic until it is recovered. See the [human source-output comparison](docs/reports/in-vivo-source-output-comparison-2026-10-11.md), [native-grid run report](docs/reports/in-vivo-base-oblique-no-regrid-2026-10-12.md), and [intermediate geometry audit](docs/reports/in-vivo-source-intermediate-geometry-audit-2026-10-12.md). PLUS integration is restartable on a committed synthetic case through iteration 26. The latest selected Java comparison (iteration 25→26) matches its tested methods exactly; the loop remains nonterminal, and the surface audit fails GM support, genus, and thickness-tail bounds. Iterations 24 and 26 return to a similar state, confirming the near two-cycle. See the [iteration-26 source and ribbon audit](docs/reports/p41z-iteration26-checkpoint-audit-2026-10-11.md) and [iteration-25 audit](docs/reports/p41z-iteration25-checkpoint-audit-2026-10-11.md). Real-human PLUS acceptance and whole-application Java equivalence remain open.

The latest source-faithful run on the labeled 65³ folded synthetic PLUS fixture completed all nine workflow stages under a 16 GB watchdog, with a 5.43 GB peak, and produced a nonterminal iteration-2 checkpoint. All 51 manifest output hashes and all 46 RAS NIfTI orientations verify. Two diagnostic NIfTIs contain source-faithful nonfinite energy values; their exact counts match the stage summaries, while membership maps remain finite. Ribbon acceptance still fails: WM/GM/CSF label support is 90.7%/71.8%/100%, meshes have 24 components and Euler characteristic 32, and surface thickness tails are out of bounds. This synthetic run confirms end-to-end execution and NaN continuation, not useful anatomy, Java equivalence, or human PLUS validation. See the [full-run report](docs/reports/plus-labeled-folded-full-run-2026-10-12.md).

Current blockers also include human source-output discrepancies and synthetic ribbon/convergence failures; these cannot be resolved by relaxing tolerances. No visualization port or FreeSurfer dependency is included. Detailed standalone-tool and workflow status is in [plan.md](plan.md); reproducible historical reports are in [docs/reports](docs/reports/).

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
