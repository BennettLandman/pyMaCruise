# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-11

The restartable BASE workflow completes on the supplied human T1 under the 16 GB memory ceiling. Orthogonal RAS regridding is enabled by default in the base and PLUS workflow commands; use `--no-regrid-inputs` for already aligned legacy inputs. The latest full-suite evidence is 1,559 passed, 4 skipped on commit `3bf4aeb2`; later changes have focused validation but no new full-suite run. See the [full-suite report](docs/reports/full-pytest-2026-10-11-lesion-dura-wrapper-scale.md) and [regridding behavior/evidence](docs/reports/orthogonal-ras-regridding-focused-2026-10-11.md).

Real-human BASE processing and selected Java stage comparisons are available, but source-output agreement and ribbon acceptance remain unresolved: held-fixed segmentation/surface thresholds fail, and the produced meshes contain many small components. PLUS integration is restartable on a committed synthetic case through iteration 24. The latest selected Java comparison (iteration 21→22) matches its tested methods exactly; iterations 23–24 remain nonterminal. Their surface audit still fails GM support, genus, and surface-thickness bounds. See the [iteration-22 selected-source audit](docs/reports/p41z-iteration22-checkpoint-audit-2026-10-11.md) and [iterations-23–24 continuation audit](docs/reports/p41z-iterations23-24-continuation-2026-10-11.md). Real-human PLUS acceptance and whole-application Java equivalence remain open.

Current blockers are upstream PLUS/source-run configuration and the synthetic ribbon/convergence failures; these cannot be resolved by relaxing tolerances. No visualization port or FreeSurfer dependency is included. Detailed standalone-tool and workflow status is in [plan.md](plan.md); reproducible historical reports are in [docs/reports](docs/reports/).

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
