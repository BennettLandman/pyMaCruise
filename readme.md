# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-07

The corrected six-stage synthetic base workflow completes on the `(172,216,164)` brain-resolution fixture: both TOADS passes converge, all 30 manifest NIfTIs are finite and RAS, and resume verifies/skips all stages (peak sampled RSS **8.526 GB** under 16 GB). Both full-grid TOADS passes now match selected VandyCRUISE Java through convergence: 12 and 33 energy ratios match at float32 precision, with voxel-exact native-grid classifications (Dice 1.0 for all 11 labels). On a `(166,208,156)` synthetic crop, MaACE now matches selected Java exactly across enhanced GM and both skeleton outputs; the comparison caught and fixed Java double-threshold arithmetic. Goldens and provenance are committed. These are synthetic stage comparisons, not whole-workflow equivalence. The selected Dura wrapper sequence matches Java on attached-dura/no-dura 41³ cases. Its one-pass stage on the actual synthetic brain-resolution base inputs meets the agreed segmentation and surface gates (23 differing mask voxels, Dice 0.9999896, Hausdorff 0.9375 mm). The selected CRUISE core Java comparison now passes inner and central surfaces but outer Hausdorff is 1.5624 mm, above the 1 mm gate; that localized deviation blocks acceptance while it is investigated. Focused restart tests pass **3/3**. The nine-stage raw-input plus workflow completes on the folded 64×80×64 fixture with verified resume. Folded TGDM Java comparison meets the agreed thresholds and folded ACE final enhanced-GM/skeleton outputs are voxel-exact. Latest full suite on `7ae1dee`: **1,220 passed, 0 failed, 0 skipped, 11 expected RAS warnings** (1,463.83s, 12.448 GB peak); it predates the latest base-workflow change. Next: isolate the outer TGDM outlier, then compare the complete CRUISE/refinement outputs and advance plus source comparisons and real-MRI validation. See the [brain-resolution Dura comparison](docs/reports/dura-brain-resolution-java-comparison-2026-10-07.md), [CRUISE core comparison](docs/reports/tgdm-core-brain-resolution-java-comparison-2026-10-07.md), [Dura wrapper comparison](docs/reports/dura-wrapper-41-java-comparison-2026-10-07.md), [Dura TGDM folded-cortex comparison](docs/reports/dura-tgdm-folded-cortex-41-java-comparison-2026-10-07.md), [MaACE brain-resolution comparison](docs/reports/ace-brain-resolution-java-comparison-2026-10-07.md), [base-run report](docs/reports/ace-indexed-heap-brain-resolution-base-run-2026-10-07.md), [full-suite report](docs/reports/full-pytest-2026-10-07-plus-uncrop.md), [folded TGDM report](docs/reports/tgdm-folded-cortex-java-comparison-2026-10-07.md), [folded ACE report](docs/reports/ace-folded-cortex-java-comparison-2026-10-07.md), and [detailed plan](plan.md).
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
