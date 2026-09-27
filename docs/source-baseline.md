# MaCRUISE Source Baseline Decision

**Decision date:** 2026-09-27  
**Status:** Selected for the first Python base-workflow port  
**Reference source:** `/Users/landmaba/dev/deleteme/masimatlab/trunk/xnatspiders/matlab/MaCRUISE/`  
**Subversion working-copy revision observed:** 9360  
**Repository URL:** `https://www.nitrc.org/svn/masimatlab/trunk/xnatspiders/matlab/MaCRUISE`

## Canonical base workflow

Use `run_maCRUISE_v3.m` as the top-level behavior and output contract, together with `pipeline/maCRUISEpipeline.m` and the functions it calls, at SVN revision 9360. Use the same checkout revision for all MATLAB reference files. The chosen VandyCRUISE algorithm source is the `VandyCRUISE` directory recorded in the project README; Python library replacements and ANTsPy integrations must reproduce this selected workflow rather than silently changing its behavior.

The project owner clarified that the version numbers on `run_maCRUISE` are bug-fix versions, and that numbered versions below V3 contain bugs. Therefore V3 is the required corrected reference, not an optional newer algorithm generation. The source also shows that its entry point sets `opts.topoCorrectVersion = 'V3'`, which activates an additional boundary-cleanup branch after multi-atlas/TOADS tissue combination. Preserve that observed behavior as part of the V3 reference contract; do not interpret the filename numbering as an algorithm-versioning or compatibility-option scheme. The entry point also exposes configurable T1 sequence, N3/N4 correction, and BET threshold, with BET as the default.

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

| Entry point | Effective topology behavior | Default TOADS skull stripping | Other notable behavior | Decision |
|---|---|---|---|---|
| `run_maCRUISE.m` | Does not set a version; `maCRUISEpipeline.m` therefore defaults to V1. | Does not set `SkullStripMethod`; `run_Toads_pipeline.m` defaults to SPECTRE. | Minimal wrapper defaults; N4 and BET threshold are supplied by downstream defaults. | Older unversioned entry point; not the port baseline. |
| `run_maCRUISE_v1.m` | Does not set a version; pipeline default is V1. | BET. | Adds configurable skull strip, T1 sequence, N3/N4 choice, and BET threshold. Its comment describes a V2-era update. | Known-buggy pre-V3 version per project owner; do not target as a supported behavior. |
| `run_maCRUISE_v3.m` | Explicitly selects V3 and triggers its boundary-cleanup branch. | BET. | Includes configurable options and a V3-specific MNI mask assignment. | **Canonical corrected base reference.** |

The supplied `_v1` and `_v3` files were compared directly. Their executable differences include the function name and the V3-specific `topoCorrectVersion` assignment; later in `_v3`, the MNI reference mask is also assigned. The plain entry point is structurally older and does not set the same option set. These source differences explain the behavior the Python port must reproduce; they do not change the project decision that the version numbers are bug-fix versions and that V3 is the corrected baseline. The exact files in the revision-9360 working copy—not recollection from comments—remain the implementation reference.

## Explicitly deferred or separate behaviors

- `run_maCRUISEplus*.m` and the `CRUISEplus*` pipelines are a separate MaCRUISEplus contract, tracked in P03 and P38–P41; they are not silently folded into the base runner.
- Do not expose pre-V3 numbered versions as supported behavior; the project owner reports that they contain bugs. The unversioned runner is also not the baseline, though its specific bug status has not been separately established.
- FreeView screenshot/report generation is optional presentation behavior, not a prerequisite for core segmentation. Surface coordinate conversion/export remains part of the output contract where downstream FreeSurfer use depends on it.
- Batch/cluster launchers and legacy/old pipelines are outside the first single-subject baseline unless later evidence shows a required scientific behavior is only present there.

## Evidence inspected

- `svn info` on the MaCRUISE directory: repository URL and working-copy revision 9360.
- `run_maCRUISE.m`, `run_maCRUISE_v1.m`, `run_maCRUISE_v3.m`, and direct `diff -u` of `_v1` and `_v3`.
- `pipeline/maCRUISEpipeline.m` and its V3 conditional boundary-cleanup branch.
- `pipeline/run_Toads_pipeline.m` defaults and `pipeline/run_maCruiseSurf_pipeline.m` V3 dispatch.

This decision validates only the selection of a source behavior contract. The project-owner clarification that numbered releases are bug fixes and pre-V3 versions are buggy is part of the decision record. It does not validate a Python implementation or numerical equivalence; those remain open roadmap items.
