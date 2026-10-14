# Case 100260 BASE execution and quality audit (2026-10-14)

## Result

The selected Python BASE workflow completed on case 100260 from its T1,
same-grid BET brain mask, and SLANT-TICV MaSeg. Its manifest reports all
**8/8 stages complete**. Integrity verification found **0 mismatches among 68
output hashes** and **0 failures among 51 finite RAS NIfTIs**. The run used
8.460 GB peak sampled RSS under a 16 GB watchdog. It ran from source snapshot
commit `326d96af750d6fdcf5acf00e174740bb861bd701` and set `flair_used_in_base_run`
to false.

This is successful workflow execution with failed anatomy-quality acceptance.
On the supplied MaSeg, GM-between-surfaces support is **65.70%** and
WM-inside-inner support is **95.70%**, below the 97% target for both regions.
Using the generated source-faithful refined segmentation for this diagnostic
raises support to **89.06% GM** and **98.21% WM**; GM still misses the target.
The three surfaces are watertight and have no nonmanifold edges, but each has
**70 connected components**, including 69 components under 100 vertices; the
largest component Euler characteristic is −94. Nesting checks find 0 central
voxels outside outer and 2 inner voxels outside central; 6,427 inner surface
voxels lie outside the supplied brain mask.

Reported diagnostic thickness is 2.507 mm mean from level-set subtraction
(275,413 crossing voxels) and 1.626 mm mean from symmetric closest-point
surface distance. These are internal measurements; no complete subject-
specific historical reference is available to assess the agreed Java/Python
surface-position gates. The supplied source log reports an unsuccessful
MaCRUISE run, and the old copied output directory is absent.

The second TOADS pass reached the configured 50-iteration cap without meeting
the `1e-4` threshold. The first pass stopped naturally at iteration 41 and
ratio `4.63128818e-05`. A separate held-fixed source replay on the saved
second-pass handoff matches selected Java exactly at the cap; see the
[TOADS Java comparison](in-vivo-100260-second-toads-java-comparison-2026-10-14.md).
That comparison is algorithm-level evidence, not whole-BASE equivalence.

## Reproduction and artifacts

The run manifest, audit inputs and generated NIfTIs remain under ignored
`build/invivo-100260-base-2026-10-13/`. Raw subject data are external and are
not committed. The integrity and quality reports are in ignored
`build/invivo-100260-run-support-2026-10-13/integrity.json` and
`quality-audit.json`.

With those run products available, rerun the output audit with:

```bash
PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 1 --report-interval-seconds 60 -- \
  env PYTHONPATH=src ./.venv/bin/python scripts/audit_in_vivo_base_products.py \
  build/invivo-100260-base-2026-10-13/audit-view-2026-10-14 \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/100260/sub-010260/ses-01/SLANT-TICVv1.2/post/FinalResult/sub-010260_ses-01_T1w_seg.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/100260/input/ses-01/anat/sub-010260_ses-01_acq-mp2rage_T1w_bet_mask.nii.gz \
  --refined-segmentation build/invivo-100260-base-2026-10-13/stages/maseg-refinement/attempt-01/SEG_refine.nii.gz \
  --dura-brain-mask build/invivo-100260-base-2026-10-13/stages/dura/attempt-01/brain_mask.nii.gz \
  --report build/invivo-100260-base-2026-10-13/audit-current-2026-10-14.json
```

The audit completed at 0.741 GB peak sampled RSS under the 8 GB watchdog.
Full-suite evidence for the current implementation is recorded separately in
the [2026-10-14 suite report](full-pytest-2026-10-14-toads-lps-axis-compatibility.md).
