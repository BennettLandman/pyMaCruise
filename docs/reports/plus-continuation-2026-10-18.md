# PLUS continuation evidence — 2026-10-18

> Historical dated record. Current testing readiness, native/Python atlas preparation and run instructions supersede stale open-work statements below; see [the October 27 readiness report](plus-testing-readiness-2026-10-27.md).

> **Owner correction, 2026-10-18:** PLUS targets `ibsr11objV3`; the earlier `ibsr8obj` requirement below is superseded. Native assets are present, but descriptor adaptation, Python NIfTI conversion, and matched validation remain open. The existing runtime commands still reflect the historical configuration and need updating after compatibility work. See [selection findings](plus-atlas-selection-2026-10-18.md).

## Final MaSeg refinement

The selected source is MaCRUISE SVN r9360. `run_maCRUISEplus_v1.m` calls
`run_SegmentationRefine` after surface production, and that wrapper sets all
four coefficients to 0.2 mm, enables central handling, and calls
`MAsegRefineBySurf`. The algorithm converts coefficients by the largest input
voxel spacing, changes cortical labels `>=100` inside the inner threshold to
local WM labels 44/45, clears cortical labels beyond the outer threshold, then
fills central 0/44/45 voxels from local cortical labels. The neighborhood is
5×5×5 and the source cap is 50 sweeps. Source `tabulate` creates rows through
the maximum positive label; the winning row index is the actual label.

The final Python full-workflow chain now appends the restartable
`maseg_refinement` stage after surfaces/thickness. It writes full-grid,
RAS-normalized `SEG_refine.nii.gz`, records input/output and implementation
hashes plus rule metrics, returns the output path, includes it in
`final_outputs`, and prints it in CLI runs. The parent hashes all Python package
sources; the standalone checkpoint fingerprints the refinement implementation
and image I/O helpers. Resume checks repair a corrupt/missing output and reject
changed inputs or implementation. The committed MATLAB fixture output was
compared label-for-label, including its full grid/affine. Failed grid checks
leave a child manifest with `status: failed` and the error.

`source_faithful=True` does not mean exact MATLAB failure behavior at image
borders. Python clips a neighborhood at the image boundary, while the MATLAB
index expression can throw for a pending edge voxel. Python also exits after a
no-change sweep, while MATLAB can repeat that sweep to 50; this second
difference leaves the fixed-point output unchanged. These deviations are
documented in the source-mapped helper and the stage evidence.

## Selected PLUS TGDM mask contract

The selected call chain is `run_CruisePlusSurf_pipeline.m` →
`run_TGDMinvcrop.m` → `run_TGDMinvcrop_jist.m` →
`MaTGDMcrop.java` (all reviewed source copies are SVN/CVS-pinned as described
below). The MATLAB/JIST argument list supplies initial membership, WM, GM,
ACE, thinned skeleton, GVF, and crop; it does not pass the optional
`Brain Mask Level Set`. `MaTGDMcrop.java` declares that `ParamVolume` optional
and forwards `null` to inner, central, and outer solves when unset. The full
Python PLUS route now leaves `brain_mask` unset for TGDM; the supplied binary
mask remains in selected preprocessing and crop use, and lower-level callers
can still explicitly request a TGDM bound.

The October 17 human comparison used a supplied signed bound in both programs;
it remains valid optional-bound evidence and does not validate the selected
wrapper. A new case-100260 null-bound Java/Python TGDM replay passes all five
surface gates on the same 140×178×134 held-fixed input arrays. The initial
level set's maximum voxel difference is 4.77e-7 with no sign changes; outer
surface Dice is 0.999994, symmetric mean distance is 0.000110 mm, and
Hausdorff is 0.225497 mm. Java peak sampled RSS is 1.259 GB under the 8 GB
watchdog. The exact comparison scope and inputs are recorded in the
[null-bound TGDM report](in-vivo-plus-tgdm-null-bound-java-parity-2026-10-18.md)
and its JSON. This remains held-fixed stage evidence, not an independent
complete PLUS workflow comparison.

A second replay was run directly from the October 17 held-fixed array bundle
(SHA-256 `c533ad2c…2edd39c9`) to exercise the newly added replay adapter. It
also passed all five surface gates; the outer surface had Dice 0.999991,
symmetric mean distance 0.000144 mm, and Hausdorff 0.710 mm. It used 2.646 GB
peak Java RSS under the 8 GB watchdog. Its detailed result is retained at
`build/human-plus-tgdm-null-bound-2026-10-18/comparison.json`. The replay adapter
now explicitly accepts the saved precomputed initial level set; treating that
signed field as membership would otherwise regenerate an empty seed.

## PLUS comparison tooling and source-composition worklist

`scripts/compare_external_plus_run.py` is a reusable, manifest-driven final
product comparator. It checks complete/terminal Python stage order, parent and
child manifest/output hashes, surface substage attempts, AtlasEM-stage
selection, shared T1/FLAIR/mask/MaSeg and atlas identity, active external
attempts/product hashes, RAS-mm grids, required labels, three level sets,
three meshes, volume thickness, and thickness summaries. Dice and surface
gates are evaluated independently; thickness/scalar values are reported
without an invented tolerance. Tool failure writes a report and exits 2;
numerical gate failure exits 1. The output map is explicit because the pinned
runtime's actual PLUS output tree has not been verified. A committed synthetic
self-comparison exercises this adapter; that is tool behavior evidence only.

| Work item | State from current source and evidence |
|---|---|
| Final refinement implementation | Implemented; MATLAB fixture labels, geometry, hashing, resume, corruption repair, changed-input/implementation rejection, and failed manifest are focused-tested. Boundary error behavior differs as documented. |
| TGDM selected optional mask | Implementation corrected and cache policy added to file-workflow settings/dependencies. Selected wrapper call chain independently inspected. Focused route/dependency tests pass. A matched null-bound source replay is reported separately below. |
| AtlasEM, initialization → registration → recrop → EM/statistics/prior updates → stopping → output → Dura | The Python composition and generated-atlas ten/eleven-stage integration exist. Selected Java/source tests cover many individual slices and optimizer paths; there is no complete original-source composed AtlasEM execution. This remains unfinished source validation, not just a missing-atlas issue. |
| LesionTOADS, Dura, crop/transfer, surfaces and thickness | Runnable composition and selected held-fixed components exist. The current human 100260 upstream products remain Python-generated; no independent selected full PLUS output has been returned. Complete numerical comparison is still missing. |
| FLAIR registration | Source contract is rigid Aladin with T1 as reference/fixed and FLAIR as floating/moving, resampled to the T1 grid. `reg_aladin` is unavailable here. ANTsPy records forward/inverse transforms and resamples to the fixed T1 reference with continuous-image interpolation, but optimizer/interpolation numerical parity is unverified. |
| External PLUS comparator/handoff | Comparator synthetic self-test passes. A specific run-record and product-map protocol is documented in [the PLUS walkthrough](../../testing/full-plus-walkthrough.md). The launch and output-layout assumptions remain explicitly untested. |

## Validation state and remaining blockers

Final combined focused validation passed **106 tests, 4 expected oblique-RAS warnings, 0 failures in 188.45 s**. The later watchdog run of the complete repository suite completed with **1,704 passed, 4 skipped, and 1 failed in 960.53 s** (30 warnings; peak sampled RSS 11.418 GB, below the 16 GB limit). The sole failure was a stale test expectation that called `ibsr8obj` selected, contrary to the owner correction. The assertion now verifies both the selected `ibsr11objV3` target and the legacy driver’s `ibsr8obj` reference; the affected test file passes **2 tests**. The full suite was not rerun after this test-only correction, so report the initial full-run count and targeted correction separately.

Installed distribution verification also passed: the sdist and wheel each contain all **107 atlas files**; the isolated installed-wheel probe checked **106 atlas files**, **35 RAS NIfTIs**, **2 runtime tables**, **7 descriptors**, and all **29** `ibsr11objV3` source files / **13** descriptor payloads. The probe reports the intended unresolved boundary: descriptor adaptation and checked RAS conversion remain open.

The earlier `ibsr8obj` atlas requirement is superseded by the owner-selected
`ibsr11objV3` package, whose native assets are present. Its descriptor needs
AtlasEM-compatible adaptation and its raw payloads need checked geometry and
voxel-order conversion; see [atlas selection findings](plus-atlas-selection-2026-10-18.md).
MATLAB and `reg_aladin` are unavailable. Therefore selected-atlas
initialization/human convergence, full original PLUS execution, Aladin
optimizer comparison, fresh selected AtlasEM composed-source validation, and
independent external-versus-Python acceptance remain open.
Do not call PLUS complete based on synthetic integration or held-fixed stages.
