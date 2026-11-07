# MaCRUISE Source Baseline Decision

**Decision date:** 2026-09-27  
**Status:** Selected for the first Python base-workflow port  
**Reference source:** `/Users/landmaba/dev/deleteme/masimatlab/trunk/xnatspiders/matlab/MaCRUISE/`  
**Subversion working-copy revision observed:** 9360  
**Repository URL:** `https://www.nitrc.org/svn/masimatlab/trunk/xnatspiders/matlab/MaCRUISE`

## N4 policy in the 0.1.6 tester release

The wrapper invokes the installed N4 command without an explicit mask. The
original runtime has not been identified. Version 0.1.5 forced a historical
ANTs 2.1.0 automatic Otsu policy; improved N4 correlation did not establish
improved final agreement. Version 0.1.6 defaults to the current installed
backend's unmasked policy (`all`) and exposes `--n4-estimation-mask legacy-otsu`
for deliberate historical-policy comparisons. The selected policy is included
in the run signature and cannot change during resume. Both modes retain the
full iteration schedule, acquisition-first N4 and native final delivery.
This supersedes the forced N4 policy described in the 0.1.5 compatibility note
below; its other source-derived corrections remain enabled.

## Adam v3.2 compatibility correction — October 29, 2026

For Adam's supplied `MaCRUISEv3.2.0` results, BASE now follows the corrected
`run_maCRUISE_v3_2_classern.m` entry point in the supplied sibling
`MaCRUISE_v3_1_0/` checkout. That directory contains a v3.2 entry point; its
name does not make this a request to support a known-buggy pre-V3 contract.
This adds the default `SepLeftRightBrain = 'SepLeftRight'` dispatch to
`pipeline/seperate_left_right_cortex.m`. Encode left/right cortical labels as
101/102, compute the existing MinFilter3D boundary, and exclude its full
support from all four cortical masks. On Adam's input, all four derived masks
then exactly match the supplied original outputs. Direct and PLUS mask callers
retain their existing label-only behavior unless explicitly requested.

Both supplied BASE source trees leave TGDM's optional Brain Mask Level Set
unset in `pipeline/run_TGDMinvcrop_jist.m`; the selected `MaTGDMcrop` Java
wrapper therefore passes null. BASE and BASE continuation now follow this
handoff. The prior Python handoff of Dura's signed field (or the original
binary skull-strip mask in continuation) added a constraint absent from that
wrapper. The optional bound remains available to explicitly configured
standalone calls. A held-bound Java/Python agreement test is evidence about
that configured component, not evidence that BASE should supply the bound.

The v3.2 topology wrapper additionally selects `background->object`,
`intensity` initialization, and `-inStop false`. The corrected revision-1.7
Java implementation ignores the supplied paint mask for that branch and
propagates from its expanded bounding box without intensity threshold stopping.
BASE now selects this upward branch, with RAS indexing retained at this
post-transfer stage. With Adam's original WM input held fixed, Python exactly
reproduces all 5,002,140 voxels of the original corrected topology crop and its
single component. Five independent selected-Java synthetic references also
match exactly. Dura and generic topology callers preserve their paint-mask
branch; the upstream Dura order and downstream topology indexing are separate
source boundaries.

BASE also restores mask-before-axial WINDOWED_SINC, source LPS order in Dura,
and the historical automatic Otsu N4 estimation policy. Full numerical
settings and native-grid final exports remain the defaults. The exact original
N4 runtime is unidentified, so these source contracts and isolated matches do
not establish complete numerical equivalence. See the
[Adam correction report](reports/adam-preprocessing-parity-2026-10-29.md).

The September decision below records the original source selection. The
October compatibility correction supersedes its missing hemisphere dispatch
and the previously implemented optional TGDM-bound handoff for BASE.

## Canonical base workflow

Use `run_maCRUISE_v3.m` as the top-level behavior and output contract, together with `pipeline/maCRUISEpipeline.m` and the functions it calls, at SVN revision 9360. Use the same checkout revision for all MATLAB reference files. The chosen VandyCRUISE algorithm source is the `VandyCRUISE` directory recorded in the project README; Python library replacements and ANTsPy integrations must reproduce this selected workflow rather than silently changing its behavior.

The project owner clarified that the version numbers on `run_maCRUISE` are bug-fix versions, and that numbered versions below V3 contain bugs. Therefore V3 is the required corrected reference, not an optional newer algorithm generation. The source also shows that its entry point sets `opts.topoCorrectVersion = 'V3'`, which activates an additional boundary-cleanup branch after multi-atlas/TOADS tissue combination. Preserve that observed behavior as part of the V3 reference contract; do not interpret the filename numbering as an algorithm-versioning or compatibility-option scheme. The entry point also exposes configurable T1 sequence, N3/N4 correction, and BET threshold, with BET as the default.

### Inputs and default entry-point settings

- Required base inputs: T1 image, MaSeg multi-atlas segmentation, working/output directory, configured atlas assets, and external scientific-processing tools selected by the Python profile. The MATLAB runner's FreeSurfer path is not a Python port requirement.
- `run_maCRUISE_v3` defaults `SkullStripMethod` to `BET`, `T1_sequence` to `T1_MPRAGE`, `N3orN4` to `N4`, `BetThreshold` to `0.1`, and `topoCorrectVersion` to `V3`. The skull-strip, sequence, correction, and threshold values can be supplied through `optsin`; the entry point sets topology correction to V3.
- The `run_maCRUISE_v3` top-level runner sets rigid-only registration (`rigOnly = 1`) for the FLAIR-related registration option. FLAIR is not a required input to this base entry point.
- `pipeline/maCRUISEpipeline.m` executes raw preparation, cortical separation, TOADS preprocessing, transfer to original space, crop, multi-atlas boundary generation/cleanup, atlas/TOADS tissue combination, V3-specific boundary cleanup, and the CRUISE surface pipeline.
- Despite an earlier `WMIsoLevel` default check, the pipeline then assigns `WMIsoLevel = 0.5` and `GVF_Lambda = 1`. It enables intermediate outputs, passes the cleaned boundary to ACE, sets `UseUpdatedCSFinACE = 2`, passes the raw cropped MA WM mask to topology correction, and enables refinement by MA. These effective values—not the earlier defaults—are the initial Python compatibility target.

### Required products to preserve

The top-level runner collects refined segmentation, raw/prepared inputs, cortical-separation products, TOADS products, uncropped level sets and scientific surface meshes, inverse-cropped TGDM inputs, thickness output, and MNI-space level sets/surfaces. It also creates FreeSurfer/FreeView-compatible ASC surfaces and FreeView batch/figure products. PyMaCRUISE ports the scientific volumes, surfaces, transforms, and thickness needed by the workflow; it does not port ASC/FreeView presentation outputs, visualization scripts, or rendering. It must not silently omit a scientific intermediate used downstream.

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
- FreeView screenshots, reports, batch scripts, rendering, FreeSurfer comparisons, and ASC presentation exports are outside this port. After acceptance, P47 is a discussion of modern Python visualization options rather than a direct visualization port.
- Batch/cluster launchers and legacy/old pipelines are outside the first single-subject baseline unless later evidence shows a required scientific behavior is only present there.

## Evidence inspected

- `svn info` on the MaCRUISE directory: repository URL and working-copy revision 9360.
- `run_maCRUISE.m`, `run_maCRUISE_v1.m`, `run_maCRUISE_v3.m`, and direct `diff -u` of `_v1` and `_v3`.
- `pipeline/maCRUISEpipeline.m` and its V3 conditional boundary-cleanup branch.
- `pipeline/run_Toads_pipeline.m` defaults and `pipeline/run_maCruiseSurf_pipeline.m` V3 dispatch.

This decision validates only the selection of a source behavior contract. The project-owner clarification that numbered releases are bug fixes and pre-V3 versions are buggy is part of the decision record. It does not validate a Python implementation or numerical equivalence; those remain open roadmap items.
