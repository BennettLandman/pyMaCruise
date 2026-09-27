# MaCRUISE Source Baseline Decision

**Decision date:** 2026-09-27  
**Status:** Selected for the first Python base-workflow port  
**Reference source:** `/Users/landmaba/dev/deleteme/masimatlab/trunk/xnatspiders/matlab/MaCRUISE/`  
**Subversion working-copy revision observed:** 9360  
**Repository URL:** `https://www.nitrc.org/svn/masimatlab/trunk/xnatspiders/matlab/MaCRUISE`

## Canonical base workflow

Use `run_maCRUISE_v3.m` as the top-level behavior and output contract, together with `pipeline/maCRUISEpipeline.m` and the functions it calls, at SVN revision 9360. Use the same checkout revision for all MATLAB reference files. The chosen VandyCRUISE algorithm source is the `VandyCRUISE` directory recorded in the project README; Python library replacements and ANTsPy integrations must reproduce this selected workflow rather than silently changing its behavior.

The reason for selecting V3 is explicit in the source: its entry point sets `opts.topoCorrectVersion = 'V3'`, whose additional boundary cleanup branch runs after multi-atlas/TOADS tissue combination. The entry point also exposes the V2-era BET default and configurable T1 sequence, N3/N4 correction, and BET threshold. This is the newest identified base entry point in the supplied MATLAB tree and has the V3 topology behavior that the earlier entry points do not select by default.

### Inputs and default entry-point settings

- Required base inputs: T1 image, MaSeg multi-atlas segmentation, working directory, and the source-era multi-atlas/MIPAV/FreeSurfer/FSL locations needed by the MATLAB runner.
- `run_maCRUISE_v3` defaults `SkullStripMethod` to `BET`, `T1_sequence` to `T1_MPRAGE`, `N3orN4` to `N4`, `BetThreshold` to `0.1`, and `topoCorrectVersion` to `V3`. The skull-strip, sequence, correction, and threshold values can be supplied through `optsin`; the entry point sets topology correction to V3.
- The `run_maCRUISE_v3` top-level runner sets rigid-only registration (`rigOnly = 1`) for the FLAIR-related registration option. FLAIR is not a required input to this base entry point.
- `pipeline/maCRUISEpipeline.m` executes raw preparation, cortical separation, TOADS preprocessing, transfer to original space, crop, multi-atlas boundary generation/cleanup, atlas/TOADS tissue combination, V3-specific boundary cleanup, and the CRUISE surface pipeline.
- Despite an earlier `WMIsoLevel` default check, the pipeline then assigns `WMIsoLevel = 0.5` and `GVF_Lambda = 1`. It enables intermediate outputs, passes the cleaned boundary to ACE, sets `UseUpdatedCSFinACE = 2`, passes the raw cropped MA WM mask to topology correction, and enables refinement by MA. These effective values—not the earlier defaults—are the initial Python compatibility target.

### Required products to preserve

The top-level runner collects refined segmentation, raw/prepared inputs, cortical-separation products, TOADS products, uncropped level sets and VTK surfaces, FreeSurfer/FreeView-compatible ASC surfaces, inverse-cropped TGDM inputs, thickness output, and MNI-space level sets/surfaces. It also creates a FreeView batch/figure product. The Python output contract should distinguish required scientific products from optional visualization products, but must not silently omit a scientific intermediate used by a downstream calculation.

The V3 runner supplies `average305_t1_tal_lin.nii.gz` and `average305_t1_tal_lin_mask.nii.gz` as MNI reference assets for the MNI output path. Asset availability and licensing are handled in roadmap step P05. Replacing the source NiftyReg/Aladin registration with ANTsPy is a later equivalence task; it does not change the reference behavior selected here.

## Differences among supplied base entry points

| Entry point | Topology selection | Default TOADS skull stripping | Other notable behavior | Decision |
|---|---|---|---|---|
| `run_maCRUISE.m` | Does not set a version; `maCRUISEpipeline.m` therefore defaults to V1. | Does not set `SkullStripMethod`; `run_Toads_pipeline.m` defaults to SPECTRE. | Minimal wrapper defaults; N4 and BET threshold are supplied by downstream defaults. | Not canonical; retain as a legacy behavior reference if needed. |
| `run_maCRUISE_v1.m` | Does not set a version; pipeline default is V1. | BET. | Adds configurable skull strip, T1 sequence, N3/N4 choice, and BET threshold. The file identifies itself with a V2 update comment. | Not canonical; use for tracing behavior before the explicit V3 topology update. |
| `run_maCRUISE_v3.m` | Explicitly V3; triggers V3 boundary-cleanup branch. | BET. | Same configurable options as `_v1`; adds the V3 topology setting and V3-specific MNI mask assignment. | **Canonical base reference.** |

The supplied `_v1` and `_v3` files were compared directly. Their executable differences are the function name and V3-specific comment/`topoCorrectVersion` assignment; later in `_v3`, the MNI reference mask is also assigned. The plain entry point is structurally older and does not set the same option set. The exact files in the revision-9360 working copy—not recollection from the comments—remain the source of truth.

## Explicitly deferred or separate behaviors

- `run_maCRUISEplus*.m` and the `CRUISEplus*` pipelines are a separate MaCRUISEplus contract, tracked in P03 and P38–P41; they are not silently folded into the base runner.
- V1 topology and SPECTRE skull stripping are not baseline behavior. They may become explicit compatibility options after the V3 path is established.
- FreeView screenshot/report generation is optional presentation behavior, not a prerequisite for core segmentation. Surface coordinate conversion/export remains part of the output contract where downstream FreeSurfer use depends on it.
- Batch/cluster launchers and legacy/old pipelines are outside the first single-subject baseline unless later evidence shows a required scientific behavior is only present there.

## Evidence inspected

- `svn info` on the MaCRUISE directory: repository URL and working-copy revision 9360.
- `run_maCRUISE.m`, `run_maCRUISE_v1.m`, `run_maCRUISE_v3.m`, and direct `diff -u` of `_v1` and `_v3`.
- `pipeline/maCRUISEpipeline.m` and its V3 conditional boundary-cleanup branch.
- `pipeline/run_Toads_pipeline.m` defaults and `pipeline/run_maCruiseSurf_pipeline.m` V3 dispatch.

This decision validates only the selection of a source behavior contract. It does not validate a Python implementation or numerical equivalence; those remain open roadmap items.
