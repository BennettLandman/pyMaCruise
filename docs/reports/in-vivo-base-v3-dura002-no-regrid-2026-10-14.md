# Case 100260 corrected-Dura no-regrid BASE run (2026-10-14)

## Result

The restartable Python BASE workflow completed on the supplied case-100260
T1, matched brain mask, and MaSeg using the selected MaCRUISE V3 Dura settings
(WM regularization 0.05, GM regularization 0.02). It used no-regrid geometry:
only axis permutations/flips to RAS, with a warning that oblique direction
cosines are discarded and voxel-to-world correspondence is not preserved.
All eight stages completed, recording 68 hashed outputs. A same-settings
`--resume` pass verified and skipped all eight stages. The final NIfTI products
are finite, RAS, and on the expected T1 grid. The workflow's RSS watchdog
sampled a 7.884 GB peak under its 16 GB limit.

Both TOADS passes reached their natural energy stopping criteria:

| Stage | Iterations | Terminal energy ratio | Stop |
| --- | ---: | ---: | --- |
| First TOADS | 21 | 4.86711651e-05 | Energy threshold reached |
| Second TOADS | 18 | 3.5033856e-05 | Energy threshold reached |

The second pass was independently rerun with the selected VandyCRUISE Java
`ToadsWorkflowSourceHarness` on the identical Dura-stripped T1, Dura mask,
atlas descriptor, geometry, and wrapper settings. The final labels match
exactly: **0 of 12,058,624 voxels differ; Dice is 1.0 for each of all 11
labels**. Java completed at iteration 18 and peaked at 3.480 GB sampled RSS.
This is selected-stage source-equivalence evidence for second-pass TOADS, not
equivalence of the complete MaCRUISE workflow.

The finished workflow's independent ribbon audit reports 85.21% MaSeg GM
support between inner/outer surfaces and 88.44% WM support inside the inner
surface. With source-faithful refined MaSeg these values are 89.64% and
96.32%. Meshes are closed and manifold but have 70 components each, including
69 components under 100 vertices; five inner voxels lie outside central and
134 inner voxels lie outside the supplied brain mask. These are anatomy and
mesh diagnostics, not independent port-acceptance criteria.

## Comparison with the supplied historical output

The provided `sample-macruise/output` tree contains surface, thickness,
refinement, and tissue products but no run log or machine-readable manifest.
Its source code/build/settings therefore remain unverified. The comparator
confirmed that its `Raw/T1.nii.gz` and `Raw/Seg.nii.gz` voxel arrays are exactly
the supplied inputs. Because the no-regrid Python and historical files have
different affine headers, this secondary comparison is by voxel index after
those exact input-array identity checks; source surface vertices are mapped
through the source Raw/T1 affine.

This historical comparison does not meet the agreed output tolerances:

- Historical `SEG_refine` differs at 37,292 voxels; 37 of 134 labels have Dice
  below 0.97, and the minimum Dice is 0.48. The fresh MATLAB/Python refinement
  comparison on fixed level sets remains a separate exact source test; this
  historical output was not produced by the controlled run used there.
- Inner, central, and outer level-set sign Dice are 0.97469, 0.97390, and
  0.97612, with 20,228, 33,398, and 47,978 sign-flip voxels.
- Symmetric mesh distances exceed the requested limits. Mean absolute
  distances are 0.144, 0.218, and 0.266 mm; Hausdorff distances are 4.129,
  5.834, and 7.059 mm for inner, central, and outer surfaces, respectively.
- Atlas GM/WM/brain-mask Dice are 0.98908/1.0/1.0; atlas brain-tissue-mask Dice
  is 0.99354.

This establishes a real-human workflow execution and a concrete discrepancy
against the supplied historical output. It does **not** establish full source
equivalence because the historical run's settings/build are unknown. It does
show that current human BASE output is not within the user's historical-output
segmentation/surface tolerances. Fresh, controlled execution of the original
full MaCRUISE source remains necessary to distinguish a Python mismatch from
different historical settings or inputs.

## Earliest historical intermediate differences

To help localize the final mismatch, the historical intermediate NIfTIs were
reoriented to closest RAS using axis permutations/flips only, then compared to
the current run on the shared 184×256×256 voxel grid. No interpolation was
used. The historical first/second TOADS files carry LPS orientation; they were
flipped to RAS before label comparisons.

| Historical intermediate vs current Python | Difference |
| --- | --- |
| N4 output | Same shape and RAS axis order, but values are not equal: MAE 1,094.82, RMSE 2,905.36, maximum absolute error 36,014.49 intensity units. Historical mean is 9,807.52; ANTsPy result mean is 8,716.97. |
| MultiAtlasStripping mask vs supplied T1 BET mask | Dice 0.95797; 111,925 differing voxels (1,357,775 vs 1,305,326 positive voxels). |
| First TOADS labels after LPS-to-RAS flips | 508,454 differing voxels; mean Dice 0.8071 across nonzero labels. The historical first-pass WM mask Dice is 0.92926. |
| Dura brain mask after LPS-to-RAS flips | Dice 0.96795; 90,404 differing voxels. |
| Second TOADS labels after LPS-to-RAS flips | 198,133 differing voxels; mean Dice 0.8814 across nonzero labels. |

These measurements identify substantial differences before the final CRUISE
surfaces. They do not isolate a Python defect: the historical run uses a
MultiAtlasStripping mask while this port uses the user-supplied BET mask, and
its N4 product differs from the selected ANTsPy-backed operation. The matched
Java/Python TOADS result above shows that TOADS itself agrees when the corrected
Python Dura handoff is held fixed. The intermediate comparison script and its
machine-readable output are
[`compare_human_base_intermediates.py`](../../scripts/compare_human_base_intermediates.py)
and the ignored `intermediate-comparison.json` in the run directory.

## Reproduction

The implementation tree is commit `9ee3ee92965a7b71f8107086dc51334f9ffd6935`;
the base-run manifest records implementation signature
`bc319ea284d76e4c0a3b10f1ce7a31d6afbe28e833f7f82d6767f220290daa9b`.
Human images, intermediates, and Java artifacts remain outside version
control or in ignored `build/` paths.

The Python run used:

```bash
PYTHONPATH=.:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/python -m pymacruise.cli base-run \
  --t1 /Users/landmaba/dev/mipav/sample-macruise/input/T1.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/input/T1_bet_mask.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/input/orig_target_seg.nii.gz \
  --output-dir build/invivo-base-v3-dura002-no-regrid-2026-10-14 \
  --no-regrid-inputs --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
  --toads-iterations 50 --dura-iterations 5 --dura-membership-level 0.3 \
  --crop-edge-width 5 --ace-diffusion-iterations 25 --gvf-iterations 200 \
  --tgdm-initial-iterations 6 --tgdm-inner-iterations 4 \
  --tgdm-central-iterations 7 --wm-iso-level 0.5 \
  --ace-use-updated-csf-in-ace 2 --source-faithful-refinement
```

The Java comparison uses the run's Dura-stripped T1 and mask, serializes the
mask-applied float32 volume in MIPAV x-fastest order, and compares its output
to the `toads-second/attempt-01/classification.nii.gz` product. The input
hashes, Java result hash, per-label comparison, full stage manifest, audit,
and historical-output comparison are in the ignored run directory named
above. The audit script is
[`audit_in_vivo_base_products.py`](../../scripts/audit_in_vivo_base_products.py);
the historical-output comparison is
[`compare_human_source_output.py`](../../scripts/compare_human_source_output.py).
