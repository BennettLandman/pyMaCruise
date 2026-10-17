# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

For human testing against the complete external MaCRUISE workflow, use the [full BASE walkthrough](docs/testing/full-base-walkthrough.md), including Linux prerequisites, virtual environment, a pinned external container, both full runs, comparisons, and return files. External execution is not yet verified on the authoring host. The [end-to-end comparison protocol](docs/testing/end-to-end-external.md) defines the acceptance boundary and remaining gaps; [Python BASE checks](docs/testing/python-base-checks.md) provide supplemental preflight checks.

## Current status — 2026-10-16

- **BASE:** The restartable eight-stage command completes synthetic and human case-100260 runs. The human run naturally converges at TOADS iterations 20/18; selected Java comparisons pass for several major stages, while historical full-output comparisons still show unresolved upstream differences. See the [human BASE run report](docs/reports/in-vivo-base-v3-no-mask-n4-source-profile-2026-10-14.md).
- **PLUS:** The restartable Python workflow runs on synthetic and human inputs. Selected Java comparisons validate multiple stages, but the selected MATLAB sequence has pre-LesionTOADS AtlasEM and Dura stages that are not yet composed into the full PLUS runner. The AtlasEM CLI now reaches Dura on generated three- and eight-class atlases. The required `ibsr8obj` descriptor/payload is missing, so selected-atlas behavior remains unvalidated. See the [source-stage audit](docs/reports/plus-source-stage-coverage-audit-2026-10-16.md) and [eight-class synthetic report](docs/reports/atlasem-eight-object-synthetic-2026-10-16.md).
- **Synthetic PLUS checkpoint:** The labeled folded-cortex case exercises source-selected surface evolution and explicit ribbon checks. Its current long TOADS profile reaches the selected iteration cap without convergence; held-fixed Java transitions and terminal TGDM comparisons pass their stated gates. This is a nonconvergent synthetic input profile, not evidence that the workflow cannot run. See the [iteration 90–99 report](docs/reports/p41z-iterations90-99-selected-java-2026-10-15.md).
- **Current validation:** Full pytest on code commit `53d8915b29fd8acc9a214d0e7685e33bac2dca69` passed **1,646 tests, 4 skipped, 0 failed**, with 21 expected orientation warnings in 728.57 seconds; peak sampled RSS was **14.788 GB** under a 16 GB watchdog. See the [eight-class AtlasEM validation report](docs/reports/atlasem-eight-object-synthetic-2026-10-16.md). This combines synthetic regression and does not imply whole-workflow source equivalence.
- **External references:** The supplied MATLAB/MIPAV/JIST runtime is unavailable on this host, and the selected `ibsr8obj` assets are not available in the staged source trees. The [source-reference run packet](source-reference-run/README.md) has exact input lists, drivers, settings, and return-file instructions for a matched original run.

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

The runtime stack is NiBabel, NumPy, and SciPy; pytest, build, and scikit-image are development dependencies, and `antspyx` is optional for ANTs-backed operations. The license/attribution inventory and follow-up TTO work are in [license_info.md](license_info.md). Currently used atlas assets are regular Git files; the selected PLUS AtlasEM `ibsr8obj` descriptor and payload are not available yet and are tracked as a P44b blocker.

## Inputs and outputs

The command-line input, assumption, and output contract is documented in [docs/cli-help.md](docs/cli-help.md). Coordinate and NIfTI rules are in [docs/image-coordinate-conventions.md](docs/image-coordinate-conventions.md). Numerical acceptance thresholds are in [docs/validation-tolerances.md](docs/validation-tolerances.md).
