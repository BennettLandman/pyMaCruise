# MaCRUISEplus Source Contract

**Decision date:** 2026-09-27  
**Reference source:** MaCRUISE MATLAB Subversion working copy, revision 9360  
**Canonical full-workflow entry point:** `run_maCRUISEplus_v1.m`  
**Shared workflow implementation:** `pipeline/maCRUISEpluspipeline.m`

## Supported MaCRUISEplus contract

The first supported plus workflow is the full MaSeg-informed MaCRUISEplus run: T1-weighted MRI + MaSeg segmentation + FLAIR image. Use `run_maCRUISEplus_v1.m` as its orchestration/output contract and `pipeline/maCRUISEpluspipeline.m` as the algorithm-stage order. The top-level runner registers FLAIR to T1 with NiftyReg Aladin in rigid-only mode, runs the shared pipeline, refines the input MaSeg from the generated inner/outer surfaces, and assembles the output products.

For the initial Python behavior profile, use the `_v1` wrapper defaults:

- `SkullStripMethod = BET`, configurable as `BET` or `SPECTRE` in the MATLAB wrapper.
- `N3orN4 = N4` is set inside `maCRUISEpluspipeline.m`.
- `LesionToadsOutput = 'dura removal inputs'` and `LesionToadsCorrect = 'false'` are set before lesion TOADS.
- `rigOnly = 1`; the MATLAB wrapper calls `reg_aladin -rigOnly` with T1 as reference and FLAIR as floating/moving image. ANTsPy is the selected Python registration tool and must reproduce the rigid transform direction, reference grid, and interpolation behavior.
- Plus CRUISE settings include `WMIsoLevel = 0.5`, `GVF_Lambda = 0.2`, intermediate output enabled, `UseUpdatedCSFinACE = 3`, `update_Enhance = 0`, and `refineTopoCorrectByMa = 1` as set by the shared plus pipeline.
- `_v1` supplies the MNI reference mask `average305_t1_tal_lin_mask.nii.gz`; its runner also carries FSL and PDF-generation configuration. FreeView/PDF figure rendering is outside the Python port scope.

`run_maCRUISEplus_v1.m` then uses the surface outputs to refine the original MaSeg (`SEG_refine.nii.gz`) and collects raw/prepared data, cortical separation, lesion-TOADS intermediates, scientific surfaces and level sets, inverse-cropped TGDM inputs, thickness, and MNI-space outputs. ASC exports, FreeSurfer/FreeView compatibility products, and figure rendering are not part of the Python port. Retain intermediate scientific files needed for reproducibility and downstream stages.

## Runner differences and interpretation

Both `run_maCRUISEplus.m` and `run_maCRUISEplus_v1.m` call the same `maCRUISEpluspipeline.m`, `run_SegmentationRefine.m`, and output-assembly sequence. The `_v1` wrapper adds an `optsin` argument, FSL configuration, PDF command selection, defaults skull stripping to BET, and supplies the MNI mask. The unversioned wrapper has no `SkullStripMethod` default, so `run_LesionToads_pipeline.m` falls back to SPECTRE; it also does not pass the MNI mask. Those differences can change outputs and are not silently merged.

The project-owner correction about lower-than-V3 bug-fix versions concerned `run_maCRUISE`; it does not by itself establish the bug status of `run_maCRUISEplus_v1`. The Python plan selects `_v1` because it is the explicitly configurable full plus wrapper and best aligns with the selected base runner’s BET/MNI-mask configuration, not because the `v1` label is presumed to mean a newer or corrected algorithm. If equivalent historical plus runs used the unversioned/SPECTRE profile, that profile can be added later as an explicit option after its need is confirmed.

## Deferred related pipelines

These are separate contracts, not aliases for the selected full plus workflow:

| MATLAB pipeline | Inputs | Distinguishing behavior | Scope decision |
|---|---|---|---|
| `maCRUISEplusOriginpipeline.m` | T1 + MaSeg + FLAIR | Uses `run_LesionToadsOrigin_pipeline`; sets `N4`, `dura removal inputs`, correction false, WMIso 0.6, GVF 0.2; omits the standard path’s `run_CleanBoundaryByToads`, and passes T1 prepared file to the surface pipeline. | Defer as an alternative backend/profile until required by users or reference runs. |
| `CRUISEplusOriginpipeline.m` | T1 + FLAIR; no MaSeg | Uses origin-space TOADS and CRUISE-surface helpers; lesion TOADS correction true; defaults to N3 in the origin helper; produces a surface workflow without the full MaSeg refinement/output assembler. | Defer; not equivalent to the selected MaSeg-informed plus run. |
| `CRUISEplusN4_pipeline.m` | T1 + FLAIR; no MaSeg | Uses standard lesion TOADS with correction false and explicit N4, then calls origin-surface processing. | Defer; distinct TOADS/surface-only entry point. |
| `CRUISEpluspipeline.m` | T1 + MaSeg + FLAIR | Direct lower-level pipeline; correction true, output `cruise inputs`, WMIso 0.6; does not itself perform full top-level registration/output assembly. | Defer; its stage choices differ from the selected full runner. |

Cluster wrappers only submit these same lower-level pipeline functions and do not define additional image-processing behavior. The standard full plus wrapper remains the only MaCRUISEplus workflow in the first implementation target. A later addition of any deferred path must get its own explicit configuration, synthetic fixtures, validation row, and reference comparison.

## Source evidence reviewed

- `run_maCRUISEplus.m` and `run_maCRUISEplus_v1.m`, including direct `diff -u` comparison.
- `pipeline/maCRUISEpluspipeline.m`, `run_LesionToads_pipeline.m`, `run_CruisePlusSurf_pipeline.m`, and `run_SegmentationRefine.m`.
- `pipeline/FlairRegistration.m` and `pipeline/aladin_reg.m`, confirming `reg_aladin -rigOnly`, reference/floating input order, and output image/transform behavior.
- `pipeline/maCRUISEplusOriginpipeline.m`, `CRUISEplusOriginpipeline.m`, `CRUISEplusN4_pipeline.m`, and `CRUISEpluspipeline.m` and their direct call sites.

This decision records the supported source behavior only. No Python plus implementation or numerical equivalence has been validated.
