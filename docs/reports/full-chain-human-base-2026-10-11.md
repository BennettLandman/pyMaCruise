# Full-chain human MaCRUISE base run — 2026-10-11

> **Superseded for V3 ACE-boundary behavior.** A later source audit found
> that this run selected ACE mode 2 without supplying the cleaned MaSeg CSF
> boundary, so the mode-2 merge was inactive. Treat the support, ACE, and mesh
> metrics below as historical diagnostics, not accepted corrected V3 output.
> The corrected tail run and current metrics are in the [boundary-enabled
> human report](in-vivo-base-tail-with-v3-boundary-2026-10-11.md).
>
> **Superseded for Dura-derived stages.** The selected MaCRUISE V3 wrapper
> passes WM/GM topology regularization 0.05/0.02. This run used Python's old
> 0.05/0.05 pair; refresh Dura, second TOADS, and downstream stages before
> treating these results as current V3 behavior. See the
> [corrected human Dura comparison](human-dura-matched-java-comparison-2026-10-14.md).

## Result and interpretation

The restartable `base-run` chain completed all eight stages on the supplied
human T1, brain mask, and MaSeg after resampling those inputs to a shared,
orthogonal RAS grid. The original oblique/sheared input grid was rejected by
the current TOADS registration implementation; the RAS regrid is recorded in
`build/invivo-base-source-v3-human-inputs-2026-10-11/ras-grid/regrid_manifest.json`.
All data and outputs remain under ignored `build/`.

This was an end-to-end human integration run using the selected V3 exposed
settings: ACE CSF mode 2 and 25 diffusion iterations, WM iso 0.5, GVF 200,
Dura five iterations, TOADS registration levels 4/50 with post-refresh 5 and
registration cap 5, TOADS cap 50, and nested TGDM initial/inner/central
counts 6/4/7. The MATLAB wrapper maps central count 7 to outer count 7; the
legacy top-level `tgdm_iterations=1` recorded in `run.json` is superseded by
those explicit values. Both TOADS passes naturally converged at 24 iterations
(energy ratios `5.83754409e-05` and `3.02124408e-05`). The run uses the
previously selected native-T1 RAS working-space path and caller-supplied brain
mask, rather than historical BET/SPECTRE preprocessing. This confirms an
end-to-end human run with the chosen stage settings, not historical
preprocessing or whole-workflow Java/MATLAB equivalence. Its Dura-derived
stages are superseded because Python used 0.05 for both topology calls instead
of the selected wrapper's 0.05 WM and 0.02 GM values.

The run manifest is `build/invivo-base-source-v3-human-rasgrid-2026-10-11/run.json`.
It records all eight stages complete, 24 final outputs, RAS geometry, and no
whole-workflow source-equivalence claim. The RSS watchdog measured a peak of
**7.470 GiB** under a **12 GiB** cap. A separate audit verified 48 stage NIfTI
hashes with no mismatches; all 66 NIfTI products inspected for orientation
were RAS. A same-settings `--resume` then verified and skipped all eight
stages; peak sampled RSS was **0.214 GiB** under a **4 GiB** cap.

## Human ribbon-quality result

The result is runnable and complete as a workflow, but it **does not pass the
agreed anatomical support target**. Against the supplied MaSeg, 528,870 of
642,291 cortical GM voxels lie between the terminal inner and outer surfaces
(**82.3412%**), and 378,979 of 428,303 WM voxels lie inside the inner surface
(**88.4839%**). These are support/containment fractions, not Dice. The
historical final segmentation needed for the ≥0.97 Dice comparison is not
available.

Source-faithful refinement leaves those absolute supported voxel counts
unchanged; its larger reported fractions (88.7522% GM and 96.4814% WM) result
from smaller label denominators. Nesting checks found 2 inner voxels outside
the central surface, 304 inner voxels outside the supplied brain mask, and no
central voxels outside the outer surface. All three extracted meshes are
edge-watertight, but each has 143 connected components, 142 of them below 100
vertices. The surface audit's mean level-set thickness is 4.1715 mm; mean
closest-point mesh distance is 2.5766 mm. Without a historical output these
thickness figures describe this run rather than demonstrate accuracy.

## Source-comparison evidence

The held-fixed Java `AlgorithmTopologyCorrection` replay on this run's exact
3,195,483-voxel fused WM-fill inputs matches Python exactly: zero differing
scalar voxels, and Dice 1.0 at all five reported thresholds. This isolates
that component only; it does not compare upstream TOADS, fusion, or the
complete MaCRUISE workflow. See
[`java-topology-comparison/report.json`](../../build/invivo-base-source-v3-human-rasgrid-2026-10-11/java-topology-comparison/report.json)
in the local ignored build output.

The direct selected Java `IsoSurfaceOnGrid` replay matches the saved Python
meshes on all three fixed human level sets: vertex/face/component counts are
identical (inner 234,222/468,444/143; central 238,058/476,116/143; outer
308,110/616,220/143). Symmetric vertex-to-triangle mean distances are
0.00000013/0.00000014/0.00000014 mm; sampled Hausdorff estimates are
0.000345/0.000488/0.000488 mm. This passes the requested 0.01 mm mean and 1 mm
Hausdorff limits for extraction on fixed fields; it does not establish anatomy
or surface-evolution equivalence. Details: [human Java surface comparison](in-vivo-full-human-jist-surface-extraction-2026-10-11.md).

The held-fixed topology-correction replay matches Python exactly on the fused
WM-fill map (0/3,195,483 differing voxels). Existing held-fixed comparisons for
fusion and TGDM also validate named components; none compares the whole
workflow.

The stage-boundary audit finds that MaSeg-label high-GM membership decreases
from 92.03% at fused Cruise input to 80.56% after ACE; terminal geometric GM
containment is 82.34%. Since thresholded map support and geometric ribbon
containment are different metrics and TGDM sub-evolution fields were not saved,
this is a localization clue rather than proof of the causal stage. See the
[stage-containment audit](in-vivo-human-stage-containment-audit-2026-10-11.md).

## Source-parameter audit

The user-selected native-RAS workflow intentionally omits BET/SPECTRE
registration and transfer. The independent source-parameter audit found the
selected high-level TOADS and Dura command values recorded, but the run manifest
does not snapshot all effective defaults. TOADS atlas coefficient 2.0,
smoothing 0.2, rigid alignment, and (18,6) connectivity are hard-coded in the
CLI rather than present in the run options. The convergence threshold and
maximum-distance values are both 1e-4, but semantic equivalence is not yet
demonstrated. The selected MATLAB Dura wrapper passes `GM2=0.02`; the Python
adapter used fixed topology regularization 0.05 for both topology passes in
this run. The older VandyCRUISE class hard-codes 0.05 for both, while the
matching MaCRUISE-selected `toads-cruise` class applies the V3 wrapper's
distinct 0.02 GM value. That source selection is now resolved and implemented;
see the [correction report](human-dura-matched-java-comparison-2026-10-14.md).
The run's `registration_max_iterations=5` is an implementation cap; source
`lastIter=0` is not separately represented in the manifest.

## Reproduction

First prepare the inputs on the orthogonal RAS target grid, as recorded in the
local regrid manifest. Then run from the repository root:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 12 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  .venv/bin/python -m pymacruise.cli base-run \
  --t1 build/invivo-base-source-v3-human-inputs-2026-10-11/ras-grid/T1_raw_ras.nii.gz \
  --brain-mask build/invivo-base-source-v3-human-inputs-2026-10-11/ras-grid/T1_bet_mask_ras.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --maseg build/invivo-base-source-v3-human-inputs-2026-10-11/ras-grid/orig_target_seg_ras_uint16.nii.gz \
  --output-dir build/invivo-base-source-v3-human-rasgrid-2026-10-11 \
  --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
  --toads-iterations 50 --dura-iterations 5 --dura-membership-level 0.3 \
  --crop-edge-width 5 --ace-diffusion-iterations 25 --gvf-iterations 200 \
  --tgdm-initial-iterations 6 --tgdm-inner-iterations 4 \
  --tgdm-central-iterations 7 --wm-iso-level 0.5 \
  --ace-use-updated-csf-in-ace 2 --source-faithful-refinement
```

Focused workflow validation after the source-profile option changes:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/pytest -q --tb=short \
  tests/test_base_resume_workflow.py \
  tests/test_base_workflow_options.py \
  tests/test_toads_cruise_core.py \
  tests/test_base_workflow_restartable.py
```

Result: **12 passed in 29.48 s**. No full pytest suite was run for this cycle;
no historical full-pipeline Java/MATLAB human output was available.
