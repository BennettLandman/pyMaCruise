# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-10

Synthetic base and plus workflows are restartable and emit RAS NIfTI products. The committed 65³ labeled folded-cortex TGDM case uses source-selected full evolution settings and passes explicit ribbon-quality checks plus all Java surface gates (Dice 1.0; maximum Hausdorff 0.070323 mm). In a regenerated held-fixed 73×91×63 TGDM stage, source-ordered JIST fast-marching matches Java across all five surfaces: zero sign differences, Dice 1.0, and maximum Hausdorff 0.0000221 mm. The natural full plus integration run was interrupted at registered-loop iteration 7 when the in-vivo base run took priority; its progress checkpoint is preserved, with peak sampled RSS 7.126 GB under a 16 GB watchdog. The in-vivo base workflow completed all seven stages on the supplied inputs, after resampling oblique RAS data onto a cardinal 1 mm RAS grid (linear T1; nearest-neighbor mask/labels) and lossless integer storage conversion for MaSeg. Both TOADS passes naturally converged in 22 and 27 iterations. The bounded smoke profile had poor tissue support (GM between surfaces 28.4%, WM inside inner 80.9%); source-selected ACE/TGDM settings improved this to 81.8% and 83.2%, with five nesting voxels, but the three meshes still each have 217 components (the largest holds 99.7% of mesh area; largest-component Euler characteristic −4). A same-settings resume verified and skipped all seven workflow stages. This is in-vivo integration evidence, not Java-equivalence, clinical validation, or a P08 pass: no reference output was available for Dice or surface-distance comparison. A spatial audit found 1,872 of 1,876 vertices on the inner mesh's non-dominant components inside the supplied mask, mostly at nearest MaSeg label 209. The mesh shells and topology need investigation. The watchdog reports the runner process, not aggregate process-tree RSS; system memory free stayed around 61–67%. See the [in-vivo report](docs/reports/in-vivo-base-workflow-smoke-2026-10-10.md), ignored `build/invivo-base-2026-10-10-ras-grid/`, and [plan.md](plan.md).

The source-selected in-vivo continuation completed five Dura iterations, a naturally converged second TOADS pass at iteration 40 (terminal energy ratio `7.55981964e-05`), the ACE/GVF/TGDM core, source-faithful MaSeg refinement, and three-surface/thickness extraction. Its 13 NIfTIs are finite and RAS; GM/WM support was 79.0%/82.7%, and each mesh has 207 components, so the case fails the current ribbon-quality criteria. Diagnostics show the five-iteration Dura mask has Dice 0.9355 with the supplied mask, with lower second-pass TOADS GM memberships at several low-support MaSeg labels; the source-correct result is unresolved. The full-resolution core peaked at 2.827 GB sampled runner RSS under a 12 GB watchdog. No historical Singularity outputs were available for numerical source comparison; this is integration evidence, not source equivalence or clinical validation. The reproducible audit and findings are in the [in-vivo report](docs/reports/in-vivo-base-workflow-smoke-2026-10-10.md) and ignored `build/`; P43 remains open.

Latest full suite: **1,505 passed, 0 failed, 0 skipped, 11 expected RAS warnings** in 1,831.87 s, with 12.464 GB peak RSS under a 16 GB watchdog, on tested implementation commit `eb04598`. Later commits `7111e4a`, `05b9f79`, `55e0bc5`, `5e94278`, and `73922eb` have direct bounded-run, Java-comparison, focused-suite, or manifest-audit evidence but are not part of that full-suite run; `c42a1a6` updates documentation only. The coarse selected TOADS source-method estimate remains **3,511–4,911 of 13,651 lines (26–36%)**; it is not a workflow-completion metric. Remaining blockers include whole-workflow source parity, broader TOADS/LesionTOADS behavior, standalone CLIs, paper annotations, and in-vivo reference/quality review. See [plan.md](plan.md), [folded-cortex validation](docs/reports/labeled-folded-cortex-validation-2026-10-09.md), and [full-suite evidence](docs/reports/full-pytest-2026-10-09-labeled-folded-cortex.md).

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
