# Same-scan human MaCRUISE source-output comparison

## Result

This is a held-fixed comparison of the downloaded Singularity outputs against saved PyMaCRUISE products on the same human scan. It is direct real-data evidence, **not** a full workflow parity or acceptance pass. The source `Raw/T1.nii.gz` and `Raw/Seg.nii.gz` arrays are voxel-identical to the supplied T1 and `orig_target_seg` inputs, respectively; all are on the same 184×256×256 oblique RAS grid.

The current Python outputs do not meet the agreed position and segmentation tolerances:

| Output | Comparison | Result | Acceptance |
|---|---|---:|---|
| Inner level set | Nonpositive-mask Dice | 0.9429 | Below 0.97 |
| Central level set | Nonpositive-mask Dice | 0.9570 | Below 0.97 |
| Outer level set | Nonpositive-mask Dice | 0.9536 | Below 0.97 |
| Inner surface | Symmetric vertex-to-triangle mean / sampled Hausdorff | 0.3058 / 6.3961 mm | Fails ≤0.01 / ≤1 mm |
| Central surface | Symmetric vertex-to-triangle mean / sampled Hausdorff | 0.3570 / 5.4603 mm | Fails ≤0.01 / ≤1 mm |
| Outer surface | Symmetric vertex-to-triangle mean / sampled Hausdorff | 0.5106 / 6.8939 mm | Fails ≤0.01 / ≤1 mm |
| Refined segmentation, `source_faithful=true` | Labels below Dice 0.97 | 95 of 139 union labels | Fails |
| Refined segmentation, `source_faithful=false` | Labels below Dice 0.97 | 92 of 134 union labels | Fails |
| Inner-surface thickness map | Positive-support Dice; MAE within shared positive support | 0.9929; 0.6054 mm | Diagnostic only; not the agreed surface-position metric |

The historical `source_faithful=true` path differs from the downloaded segmentation at 107,094 voxels (median/mean/minimum label Dice **0.9520/0.9103/0.0**); labels 1, 2, 3, 5, and 6 appear only in that Python result. This is consistent with the known MATLAB `maxind` row-index bug. With the corrected `source_faithful=false` path, those five extra labels disappear and no labels are unique to either image, but 80,939 voxels still differ (median/mean/minimum label Dice **0.9559/0.9484/0.5185**) and 92 of 134 union labels remain below 0.97. Thus the intentional correction explains the label-set mismatch, but not the broader source-refinement mismatch. Both modes are reported because the source run's exact compiled version/settings were not captured.

The source and Python level-set arrays were compared by voxel index. The source level-set NIfTI headers do not retain the raw T1 oblique affine, although the arrays have the same shape. The source VTK points are voxel indices; they were transformed through the source Raw/T1 affine into RAS millimeters before surface metrics. Surface metrics therefore compare the meshes in the same physical frame. Hausdorff is the existing sampled-vertex metric, not a continuous exact Hausdorff computation.

## Source tissue-mask provenance

The cached source `Atlas_GM` contains 628,552 positive voxels, while rebuilding its mask from the downloaded `Raw/Seg` and the source V3 label selection gives 642,437. Exactly 13,885 expected GM voxels are absent, with no extra positive voxels. `Atlas_BrainTissueMask` has the same missing set. The downloaded WM mask and `Atlas_Mask` each match their source-label reconstruction exactly (428,225 voxels).

The mismatch is not explained by orientation, scaling, or the supplied BET mask. Source preparation code can skip existing outputs without checking input hashes, and no source run log or machine-readable settings/input manifest was present in the downloaded tree. Treat the GM and brain-tissue mask files as provenance-uncertain; do not use them as authoritative parity references until regenerated or their provenance is established. The level-set and final segmentation differences above are still useful comparison evidence, but exact source run settings remain unverified.

The Python tissue masks compared by the script reconstruct to 642,437 GM and 428,225 WM voxels. Against the downloaded cached masks, binary Dice is 0.9891 for GM, 0.9935 for `Atlas_BrainTissueMask`, and 1.0 for both WM and `Atlas_Mask`.

## Reproducibility

The comparison script is [`scripts/compare_human_source_output.py`](../../scripts/compare_human_source_output.py). Run it from the repository root with the downloaded source output and input directories available at their current locations:

```bash
PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 \
  --interval-seconds 1 --report-interval-seconds 60 -- \
  .venv/bin/python scripts/compare_human_source_output.py
```

The script writes detailed per-label, voxel, topology, and thickness metrics to the ignored `build/source-output-direct-comparison-2026-10-11/comparison.json`. It does not copy human imaging data into the repository. Python level-set, mesh, and thickness comparisons use saved exact-source-grid products at the defaults in the script. The source-faithful segmentation numbers are from the script's saved result. The corrected refinement run is a second invocation of the same command with `--no-source-faithful`, writing `SEG_refine_corrected.nii.gz` and `SEG_refine_corrected-summary.json` in the same ignored build directory; its label metrics were computed against the downloaded source `SegRefine/SEG_refine.nii.gz`. These are held-fixed stage comparisons rather than a fresh end-to-end rerun.

Run completed with exit status 0 under the 8 GB watchdog (peak sampled RSS **1.692 GB**); the script's zero exit means the comparison was produced, not that the acceptance gates passed. `surface_gates_pass` is false.

## Interpretation and next action

The large surface differences are consistent with the substantial level-set sign disagreements; they are not a small mesh-extraction-only discrepancy. Next, trace the level-set and refinement differences using the source run's exact settings and regenerate source masks from the verified Raw/Seg pair if the original pipeline can be rerun. Preserve the current mismatch as a baseline; do not change tolerances to make it pass. A human source-output pass requires every required label to reach Dice ≥0.97, symmetric mean surface distance ≤0.01 mm, and Hausdorff distance ≤1 mm, with source settings and image-space mapping accounted for.
