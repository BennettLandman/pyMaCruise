# Human TOADS membership comparison by voxel index — 2026-10-12

## Scope

This compares the source Singularity output's first- and second-pass TOADS
membership maps with the selected-profile Python BASE run. The source maps are
`(184, 256, 256)` LPS images; the Python maps are on the orientation-only RAS
grid with the same shape. For array comparison, source axes 0 and 1 were
flipped to align LPS voxel order with RAS. Source T1 and MaSeg voxel arrays
were independently verified identical to the provided input arrays.

The saved source output includes `ToadsSeperate/MultiAtlasStripping`, although
`run_maCRUISE_v3.m` defaults `SkullStripMethod` to BET unless the caller
overrides it. Its saved `MAstrip_mask` is exactly the binarized input MaSeg
(Dice 1.0, zero differing voxels). The user-supplied `T1_bet_mask` differs
from that source mask on the shared raw-T1 grid (Dice **0.95797**, 111,925
voxels differ). `run_MultiAtlasStripping.m` then sends the masked T1 through
`run_Spectre_jist.m`, whose SPECTRE2010 invocation requests 12-DOF affine
registration, trilinear interpolation, and `inResample=false`. The source run
therefore contains both a different brain-mask construction and a distinct
pre-TOADS spatial operation. The output folder establishes that the
MultiAtlas stage ran, but no saved run options/log establish how the wrapper
overrode its BET default. The Python run intentionally remains in the supplied
T1 working space and uses the user mask, per the selected workflow design.

The source TOADS NIfTI headers use an LPS diagonal and zero origin, while
Python outputs use an axis-aligned RAS affine anchored at the input center.
Consequently these are **voxel-index membership diagnostics**, not physical
grid equivalence. They do not isolate how much of the difference is
attributable to the source MaSeg-derived mask, its SPECTRE preprocessing,
registration settings, or subsequent TOADS updates. These findings do not
justify changing the selected native-space Python workflow or substituting
MaSeg for the user-supplied brain-mask input.

## Measurements

| Pass/map | Mean absolute membership difference (0–255 scale) | Sampled correlation | Dice at membership > 102/255 | Source / Python voxels above threshold |
|---|---:|---:|---:|---:|
| First GM | 4.0971 | 0.9350 | 0.90268 | 868,204 / 798,573 |
| First WM-fill | 1.7313 | 0.9779 | 0.93328 | 512,635 / 455,624 |
| First WM-mask | 1.2949 | 0.9279 | 0.92926 | 455,086 / 410,577 |
| Second GM | 3.3200 | 0.9401 | 0.91094 | 869,880 / 779,528 |
| Second WM-fill | 1.1803 | 0.9876 | 0.95587 | 484,680 / 452,503 |
| Second WM-mask | 0.8404 | 0.9511 | 0.95257 | 430,061 / 407,907 |

The selected-profile level-set Dice reaches 0.97102–0.97240, but the membership
maps already differ materially before surface evolution, particularly in
thresholded cortical GM support. This is evidence of an upstream human
workflow gap, not proof of which upstream operation causes it. The JIST
surface extractor remains matched on fixed Python level sets; the residual
surface-position error is consistent with upstream field differences.

## Reproduction

The comparisons use source maps under
`/Users/landmaba/dev/mipav/sample-macruise/output/ToadsSeperate/` and the
selected-profile run under ignored `build/`:

```text
First-pass Python maps:
  build/invivo-base-orientation-only-source-profile-2026-10-12/stages/toads-first/attempt-01/dura-inputs/t1-n4_{gm,wmfill,wmmask}.nii.gz
Second-pass Python maps:
  build/invivo-base-orientation-only-source-profile-2026-10-12/stages/toads-second/attempt-01/dura-inputs/stripped_t1_{gm,wmfill,wmmask}.nii.gz
Source maps:
  sample-macruise/output/ToadsSeperate/Toads1stPass/*toads1{gm,wmfill,wmmask}.nii.gz
  sample-macruise/output/ToadsSeperate/Toads2ndPass/*toads2{gm,wmfill,wmmask}.nii.gz
```

Reproduce the mask and membership measurements with:

```sh
PYTHONPATH=src .venv/bin/python scripts/compare_human_toads_membership.py \
  --python-run build/invivo-base-orientation-only-source-profile-2026-10-12 \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/input/T1_bet_mask.nii.gz \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/input/orig_target_seg.nii.gz \
  --report build/human-toads-membership-source-profile-2026-10-12/comparison.json
```

The Python and source affine headers are intentionally not treated as a shared
world coordinate system. No changes were made to source files or image data.
