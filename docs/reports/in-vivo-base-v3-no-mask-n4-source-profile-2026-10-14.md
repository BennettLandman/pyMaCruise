# Case 100260 BASE with the corrected unmasked N4 profile — 2026-10-14

## Execution and integration result

The eight-stage no-regrid BASE workflow completed from the supplied T1, its
matching brain mask, MaSeg, and bundled cruise atlas descriptor. It used the
V3 Dura topology pair (WM `0.05`, GM `0.02`) and the corrected unmasked ANTs
N4 profile (`rescale_intensities=False`). Both TOADS passes stopped
naturally: first pass at iteration 20 (`4.29838547e-5`) and second pass at
iteration 18 (`4.44517791e-5`). This is a real-human integration run; it does
not by itself establish whole-workflow source equivalence or anatomical
acceptance.

The run manifest records all **8/8 stages complete**, **68/68 output hashes
matching**, and **51/51 finite RAS NIfTIs**. A same-settings `--resume` call
verified and skipped all eight stages. The full-run RSS watchdog peaked at
**8.027 GB** under a 16 GB limit; the resume check peaked at 0.700 GB under
an 8 GB limit. Implementation was commit
`d4b08f8a023a869fa117a7cc1afdfcd317c8a82c`. The run's
`implementation_signature` is
`dfc95ed2578c784c771cecbbcd0f1ce55c629edfd0db416b2e1065e80dded9b6`.

## Selected Java comparison for both TOADS handoffs

Each source replay used the exact current Python intensity volume and mask,
the same source atlas descriptor, 1 mm RAS axis grid, V3 wrapper registration
settings, and a 50-iteration maximum. Java raw labels were reshaped in
Fortran order to MIPAV's x-fastest voxel order and compared to the matching
Python classification. Both source calls stopped at the same natural
iteration as Python:

| Handoff | Python / Java stop | Changed voxels | Label Dice | Java peak RSS |
| --- | ---: | ---: | --- | ---: |
| First TOADS | 20 / 20 | 0 / 12,058,624 | 1.0 for all 11 labels | 3.407 GB |
| Second TOADS after corrected Dura | 18 / 18 | 0 / 12,058,624 | 1.0 for all 11 labels | 4.507 GB |

Machine comparisons and held-fixed Java inputs are saved under the ignored
run directory in `java-source-first/` and `java-source-second/`. Java result
hashes are `c9d8ccd8705fe7d58b59a7ebfb29f98951bf97fa26d2382e6231a08e00a8dfa3`
and `4ded764f4ac99bff7c76049c5395c9176bba588074b0d3d3da2ff00b63e6611b`.
This is selected-source equivalence for both TOADS stages on matched inputs,
not a complete MATLAB/JIST BASE replay.

## Anatomy, meshes, and thickness diagnostics

Using the original supplied MaSeg, GM-between-surfaces support is **85.23%**
and WM-inside-inner support is **88.42%**. With the Python source-faithful
refined MaSeg used for this diagnostic, support is **89.66% GM** and **96.32%
WM**. The inner/central/outer meshes are closed by edge incidence and have no
nonmanifold edges; each has **67 connected components**, including 66
components under 100 vertices. Four inner-mask voxels lie outside central,
zero central-mask voxels lie outside outer, and 132 inner voxels lie outside
the supplied brain mask. These remain diagnostics, not a requirement to
improve the original algorithm.

Level-set thickness has mean **4.226 mm** over 604,953 crossing voxels;
symmetric closest-surface thickness has mean **2.653 mm**. These are two
different measurements and neither replaces source comparison.

## Supplied historical output comparison

The downloaded historical output has no run manifest/log, and its upstream
mask/registration path differs from this owner-selected workflow. It is a
secondary diagnostic, not a controlled source oracle. The source `Raw/T1`
and `Raw/Seg` arrays still match the provided input files voxel-for-voxel;
files were compared by voxel index after the explicitly requested axis-only
RAS handling. A nearby case-100260 v3.2.0 log does not supply provenance for
this output: it names a different 176×240×256 BIDS MP2RAGE/SLANT input, records
MultiAtlas skull stripping, and ends in status 1 with `SEG_refine.nii.gz`
missing; the associated versioned `Output` tree has no files. The downloaded
root `output/Raw` arrays are 184×256×256 and exactly match the separate
provided `sample-macruise/input` files. Thus the root output remains a useful
secondary comparison, but its run settings are still unknown.

- Historical N4 vs current no-mask ANTsPy N4: whole-volume MAE **6,165.91**,
  RMSE **14,573.26**, max absolute difference **137,189.73** intensity units.
  The historical output range is `[0, 210566.03]`; current Python is
  `[0, 73851.80]`. This records a real mismatch; it does not override the
  wrapper's effective no-mask rescaling behavior or identify the historical
  executable/settings.
- Historical MultiAtlas stripping mask vs supplied BET mask: Dice **0.95797**
  with 111,925 differing voxels. Therefore downstream historical TOADS
  differences are confounded by distinct input masks and N4 products.
- Historical vs Python first and second TOADS classifications differ at
  509,194 and 198,646 voxels respectively after axis-only RAS comparison.
  On held-fixed current inputs, both TOADS stages match selected Java exactly.
- Historical refined segmentation differs at **37,372 voxels**; **35/134**
  union labels are below Dice 0.97, mean Dice is 0.97262, and minimum Dice is
  0.44.
- Historical level-set sign Dice is **0.97459 inner**, **0.97384 central**,
  and **0.97611 outer**. Surface symmetric mean distances are **0.14414,
  0.21887, 0.26637 mm** and Hausdorff distances are **4.138, 5.831,
  7.059 mm** (inner, central, outer). All three surfaces fail the agreed
  `0.01 mm` mean and `1 mm` Hausdorff gates. These values do not isolate the
  mesh extractor from level-set and upstream differences.

Thus current evidence is: complete restartable human execution; exact
selected Java parity on both TOADS handoffs; and a failed, provenance-limited
historical whole-output comparison. Current BASE equivalence remains open.

## Reproduction

The current output is ignored under
`build/invivo-base-v3-dura002-no-regrid-n4off-2026-10-14/`; no subject data,
Java arrays, or run outputs are committed. Run from the repository root:

```bash
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  MPLCONFIGDIR=build/test-mpl .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/python -m pymacruise.cli base-run \
  --t1 /Users/landmaba/dev/mipav/sample-macruise/input/T1.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/input/T1_bet_mask.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/input/orig_target_seg.nii.gz \
  --output-dir build/invivo-base-v3-dura002-no-regrid-n4off-2026-10-14 \
  --no-regrid-inputs --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
  --toads-iterations 50 --dura-iterations 5 --dura-membership-level 0.3 \
  --crop-edge-width 5 --ace-diffusion-iterations 25 --gvf-iterations 200 \
  --tgdm-initial-iterations 6 --tgdm-inner-iterations 4 \
  --tgdm-central-iterations 7 --wm-iso-level 0.5 \
  --ace-use-updated-csf-in-ace 2 --source-faithful-refinement
```

Run the same command with `--resume` to verify the recorded stage hashes and
skip completed stages. Focused N4 regression evidence is recorded in the
[N4 profile correction](n4-no-mask-rescaling-default-correction-2026-10-14.md).
No full test-suite evaluation was performed as part of this human run.
