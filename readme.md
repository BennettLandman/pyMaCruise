# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-11

The restartable human BASE run completed from hash-verified 1-mm RAS T1, user brain mask, and MaSeg. All seven stages passed; TOADS stopped naturally at iterations 24 and 26; peak sampled RSS was 7.302 GB under a 16 GB watchdog. Selected JIST transfer matches Python bitwise for all three 19.6-million-voxel maps, the selected cropper matches bounds and all seven maps exactly, and attempt-02 selected-JIST surface extraction matches the Python meshes (136 components per surface; symmetric mean distance below 0.000001 mm). The island topology is already present in the seed: its negative mask exactly matches the topology-corrected WM-fill threshold (408,915 voxels, XOR 0; 82 singleton components). A bounded TGDM-only replay matches saved inner/central/outer outputs bitwise on all 3,130,400 source-written crop voxels; singleton counts evolve from 82 to 127 in the delivered inner field. The anatomical validity of those components remains unresolved. The output audit found 82.384% GM and 88.442% WM label support (not Dice). The supplied Singularity output directory is currently empty, so the agreed segmentation/surface tolerances and whole-workflow equivalence remain unverified. See the [attempt-02 human BASE report](docs/reports/in-vivo-base-rasgrid-run-attempt-02-2026-10-11.md), [attempt-02 surface comparison](docs/reports/in-vivo-base-attempt-02-jist-surface-comparison-2026-10-11.md), [human TGDM island lineage](docs/reports/human-tgdm-island-lineage-2026-10-11.md), [ACE correction report](docs/reports/in-vivo-human-ace-ordered-fmm-correction-2026-10-11.md), and [detailed plan](plan.md).

The TOADS constructor-border and loop-carried skeleton state follows selected Java behavior. Full suite on implementation commit `9cade54`: **1,531 passed, 0 failed, 0 skipped, 11 expected RAS warnings** in 801.22 s; peak sampled RSS was 13.450 GB under the 16 GB watchdog. The selected TOADS source-line estimate remains **3,511–4,911 of 13,651 lines (26–36%)**; no full suite has been run for newer P41z changes. P41z iteration-6 thinning matches Java across all **267,638 queue events**; its 13-class statistics, membership/objective/energy ratio, rigid-registration trace, transform, 13 support bounds, wrapper stop/grow choice, and relation refresh match selected Java on held-fixed synthetic states. With the Java-validated post-thinning ordering and class volumes, the iteration-6 growth event trace, final state, and post-growth relation refresh also match the selected Java harness byte for byte (**1 passed in 116.56 s; 0.700 GB peak under a 4 GB watchdog**). The prior objective still comes from a hash-verified Python integration checkpoint, not an independent Java iteration-5 replay. The complete plus run remains nonterminal at energy ratio 0.131597 vs 0.0001 and fails ribbon-quality checks. Synthetic component matches, selected-stage source comparisons, and human BASE integration remain distinct evidence categories. See the [human BASE report](docs/reports/in-vivo-base-rasgrid-run-attempt-02-2026-10-11.md), [latest full-suite report](docs/reports/full-pytest-2026-10-11-toads-skeleton-carry.md), [P41z membership checkpoint](docs/reports/p41z-iteration6-membership-carry-2026-10-11.md), and [detailed plan](plan.md).

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
