# In-vivo base MaCRUISE smoke run — 2026-10-10

## Result

The complete seven-stage base workflow ran to `status=complete` on the supplied in-vivo T1, user-provided brain mask, and MaSeg. This was a bounded integration run, not a source-default numerical validation: CRUISE core used one TGDM iteration per evolution and one ACE diffusion iteration. The resulting surfaces fail the tissue-support quality checks, so this run demonstrates executable in-vivo integration but does **not** pass base-workflow quality acceptance.

No source images or generated image products were copied into Git. Inputs and outputs remain under the ignored `build/` tree.

## Input preparation

- No FLAIR was provided, so the base MaCRUISE route was selected.
- The user-provided `T1_bet.nii.gz` and matching `T1_bet_mask.nii.gz` were used; no BET or FreeSurfer operation was run.
- The images were RAS-oriented but oblique. The supplied TOADS geometry adapter rejects oblique/sheared grids, so local working copies were resampled to a cardinal 1 mm RAS grid. T1 used linear interpolation; mask and MaSeg used nearest-neighbor interpolation. Grid dimensions changed from `(184, 256, 256)` to `(208, 309, 305)` to cover the source voxel-cell bounds.
- `orig_target_seg.nii.gz` had float64 storage but finite integral labels. It was stored as uint16 without interpolation or label-value changes before grid resampling. All 135 label values were retained.
- The regrid inputs and transform/label audit are in ignored `build/invivo-base-2026-10-10-inputs/`; no files outside the repository were changed.

## Run settings and stage results

Run directory: `build/invivo-base-2026-10-10-ras-grid/`. The `run.json` contains stage commands, hashes, options, atlas resources, and RAS conventions.

| Stage | Result |
|---|---|
| N4 | Complete |
| First TOADS | Natural convergence in 22 iterations; terminal energy ratio `9.07685535e-05` |
| Dura | Complete; one configured evolution iteration |
| Second TOADS | Natural convergence in 27 iterations; terminal energy ratio `2.59382614e-05` |
| CRUISE core | Complete; ACE 1 iteration, GVF 200 iterations, TGDM initial/inner/central/outer 1/1/1/1 |
| MaSeg refinement | Complete; source-faithful option enabled |
| Surface and thickness | Complete |

The full run command used a four-level TOADS registration pyramid, 50 registration and TOADS iteration caps, a post-refresh registration cap of 5, Dura 1, ACE 1, GVF 200, and TGDM 1. Both TOADS passes stopped by their energy criteria before their iteration caps.

## Output audit

An independent audit verified all **46** stage-output SHA-256 entries. All **31** output NIfTIs were finite, had the manifest shape, and were RAS-oriented. The top-level manifest reports all seven stages complete, `converged=true`, and explicitly makes no full Java-equivalence claim. The 24 linked final products include the RAS level sets, refined segmentation, three RAS-mm VTK surfaces, inner-surface thickness, level-set thickness, summaries, and provenance.

The thickness summary reports mean surface closest-point distance **1.827 mm** and mean level-set subtraction thickness **2.033 mm**. All three meshes are finite, watertight by edge incidence, and have zero nonmanifold edges, but each contains **217 connected components**. These topology results are not acceptable evidence of a single clean cortical ribbon.

Using the selected MaCRUISE GM/WM MaSeg label sets and zero-level regions:

- GM voxels between inner and outer surfaces: **182,215 / 642,291 (28.4%)**.
- WM voxels inside the inner surface: **346,571 / 428,303 (80.9%)**.
- Inner voxels outside central region: **21**.
- Central voxels outside outer region: **0**.

The strong support deficits and fragmented meshes fail the current qualitative ribbon checks. The shortened TGDM settings are a leading explanation; the report does not claim they are the sole cause.

## Source-selected CRUISE core and surface follow-up

To test whether the smoke-run quality deficits came from shortened evolution, the converged TOADS/Dura products from the first run were reused and the CRUISE core was rerun with ACE diffusion 25 and source-selected TGDM iteration settings: initial 6, inner 4, central 7, and outer 7 (the selected wrapper maps central iterations to outer iterations). Source-faithful MaSeg refinement and surface/thickness then ran on those new level sets. Upstream registration caps and the one-pass Dura setting remain those of the smoke run; this is not a source-default whole-workflow run.

All **22** files from this core/refinement/surface follow-up passed SHA-256 recording in the ignored audit; all **13** NIfTIs were finite and RAS-oriented. Tissue support improved:

- GM between inner and outer surfaces: **525,643 / 642,291 (81.8%)**.
- WM inside the inner surface: **356,149 / 428,303 (83.2%)**.
- Inner voxels outside central: **5**; central outside outer: **0**.
- Mean symmetric surface thickness: **2.850 mm**; mean level-set subtraction thickness: **4.556 mm** (different methods; not a parity comparison).

All three meshes are watertight by edge incidence with no nonmanifold edges, but each still has **217 components**. The largest component contains 99.7% of total triangle area; 215 components have fewer than 100 vertices, and the largest component Euler characteristic is `−4`. Thus one dominant surface is present alongside many small shells and nontrivial topology. The support metrics improve over the smoke run, but no source/reference images were available to compute the agreed per-label Dice or surface-distance tolerances. The in-vivo case is not declared numerically equivalent or fully quality-accepted.

A follow-up spatial audit of the inner mesh found that **1,872 of 1,876 vertices** belonging to its 216 non-dominant components map inside the supplied brain mask. The nearest MaSeg voxel labels for these vertices were 209 (1,569), 11 (246), 35 (22), 0 (21), 122 (12), and 134 (6). This rules out the simple explanation that nearly all small shells are outside-mask fragments; the label meanings and whether the shells are valid structures still need investigation. This diagnostic used nearest-voxel label lookup and is not an anatomical classification.

The machine-readable audit and source-profile outputs are in ignored `build/invivo-base-2026-10-10-ras-grid/source-selected-core-audit.json` and `build/invivo-base-2026-10-10-ras-grid/source-selected-core/`.

## Source-selected full-resolution continuation

After the first source-selected ACE/TGDM follow-up, the same prepared in-vivo case was carried through the higher-iteration upstream stages and a fresh CRUISE core/tail. The first TOADS pass reused the natural 22-iteration result from the complete base smoke workflow. The second Dura evolution ran five configured iterations using the source-selected Dura parameters, and the second TOADS loop then naturally converged in **40 iterations** with terminal energy ratio **`7.55981964e-05`**; it did not reach the 50-iteration cap. Its full stage log and products are in ignored `build/invivo-base-2026-10-10-ras-grid/source-selected-upstream/`. The Dura reference-gradient and non-flipping TGDM vectorizations added for this run were checked against their scalar/source implementations on code commit `73922eb`: `PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl ./.venv/bin/pytest -q tests/test_tgdm.py tests/test_toads_dura.py` — **55 passed**, two expected RAS warnings, 25.81 seconds.

The source-selected continuation’s computational code was commit `73922eb`; commit `c42a1a6` changed only README/plan documentation while the surface job was running.

The fresh CRUISE core used ACE diffusion 25, GVF 200, and TGDM initial/inner/central iterations 6/4/7, with the selected MaSeg cerebellum removal and a two-voxel crop edge. Source-faithful MaSeg refinement completed, and the source-derived JIST connectivity-18 surfaces plus thickness products were generated. The full-resolution core took substantially longer than the earlier profile; its direct runner peak sampled RSS was **2.827 GB** under a 12 GB watchdog. The new audit inventory contains **22 hashed products**, including **13 finite RAS NIfTIs**; no invalid or non-RAS NIfTI was found.

The new profile did **not** pass the agreed ribbon-quality criteria:

- GM MaSeg label support between inner and outer surfaces: **507,392 / 642,291 (79.0%)**, below the 97% minimum.
- WM MaSeg label support inside the inner surface: **354,294 / 428,303 (82.7%)**, below the 97% minimum.
- Nested regions: **1** inner voxel outside central; **0** central voxels outside outer.
- Each surface is edge-watertight with no open or nonmanifold edges, but each has **207 connected components**, including **205 components under 100 vertices**. The largest component Euler characteristic is **−2** for all three meshes.
- Symmetric closest-surface thickness mean: **2.851 mm** (SD **0.962 mm**). Level-set subtraction thickness mean: **4.455 mm** (SD **1.269 mm**).
- MaSeg labels **208 (91,275 voxels)** and **209 (24,715 voxels)** remain outside the selected label groups through 207; their anatomy/role is not inferred.

An internal comparison against the earlier source-selected core (which reused the bounded upstream outputs) shows the fuller upstream profile did not improve tissue support: GM changed from 81.84% to **79.00%** (−2.84 percentage points), WM from 83.15% to **82.72%** (−0.43 points), while inner-outside-central nesting improved from five voxels to **one**. The largest per-label support loss was label 161 (71.4% to **44.3%**, −27.2 points), followed by labels 135 (−17.7 points), 125 (−12.2), and 122 (−11.8). This profile comparison is diagnostic only; it is not a source-output comparison and does not identify the anatomical label names.

An upstream diagnostic shows the five-iteration Dura mask has **1,185,356** voxels versus **1,305,385** in the supplied mask (Dice **0.9355**; 140,290 supplied-mask voxels were removed, 20,261 added). The second-pass TOADS GM map’s mean over selected MaSeg GM labels is **203.32** versus **211.81** in the earlier profile; the fraction above the diagnostic scaled threshold `0.4 × 255 = 102` is **86.16%** versus **89.82%**. Label 161’s fraction above that threshold changed from 75.6% to **46.5%**, label 135 from 82.4% to **63.4%**, label 125 from 77.1% to **60.8%**, and label 122 from 85.0% to **68.9%**. This localizes a substantial profile-dependent change to the Dura/TOADS handoff before CRUISE, but does not establish which output is source-correct; a matched reference remains necessary.

The 206 non-dominant components on the new inner mesh contain **1,752 vertices**; **1,748** nearest-map inside the supplied brain mask. Of those in-grid shell vertices, **1,452** map to MaSeg label 209, 252 to label 11, 15 to label 35, 6 each to labels 122 and 208, and 21 to background. The same earlier shell audit also concentrated most small components at label 209. Label 209 is not named in the historical MaCRUISE label-name table, so its anatomy and whether these islands should be included remain unresolved.

These are in-vivo integration and quality measurements, not evidence of Java/MATLAB equivalence or clinical validity. The agreed 3% Dice and 1 mm Hausdorff thresholds still require matched reference segmentations/surfaces; the 0.01 mm position criterion also requires the historical source mesh. No such outputs were present in `sample-macruise/output` at inspection.

The audit is reproducible with:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl ./.venv/bin/python scripts/audit_in_vivo_base_products.py \\
  build/invivo-base-2026-10-10-ras-grid/source-selected-full-toads-core \\
  --maseg build/invivo-base-2026-10-10-inputs/ras-grid/orig_target_seg_ras_uint16.nii.gz \\
  --brain-mask build/invivo-base-2026-10-10-inputs/ras-grid/T1_bet_mask_ras.nii.gz \\
  --watchdog-log build/invivo-base-2026-10-10-ras-grid/source-selected-upstream/logs/cruise-core-second-toads.log \\
  --comparison-core-dir build/invivo-base-2026-10-10-ras-grid/source-selected-core \\
  --toads-gm build/invivo-base-2026-10-10-ras-grid/source-selected-upstream/toads-second/dura-inputs/stripped_t1_gm.nii.gz \\
  --comparison-toads-gm build/invivo-base-2026-10-10-ras-grid/stages/toads-second/attempt-01/dura-inputs/stripped_t1_gm.nii.gz \\
  --dura-brain-mask build/invivo-base-2026-10-10-ras-grid/source-selected-upstream/dura/brain_mask.nii.gz \\
  --report build/invivo-base-2026-10-10-ras-grid/source-selected-full-toads-core-audit.json
```

Audit output and all derived images remain under ignored `build/`; no human image data were committed.

## Evidence limits and follow-up

This is one in-vivo integration experiment with a documented oblique-grid resampling step. It is not comparison against the historical Singularity output: `sample-macruise/output` contained no reference products at inspection. It is also not Java equivalence, source-default execution, clinical validation, or a test of MaCRUISEplus.

The source-selected CRUISE core/tail follow-up is complete. A same-settings resume of the original seven-stage workflow verified and skipped all seven stages. Remaining in-vivo validation needs a reference run/output for per-label Dice and symmetric surface-distance comparison, plus investigation of the extra surface shells and Euler characteristic.

The RSS watchdog used by the smoke command samples its direct workflow-runner process, not the full process tree. System-wide `memory_pressure` reported 66–69% free during/after the run; the 7.053 GB runner RSS peak must not be interpreted as total process-tree peak RSS.
