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

## Evidence limits and follow-up

This is one in-vivo integration experiment with a documented oblique-grid resampling step. It is not comparison against the historical Singularity output: `sample-macruise/output` contained no reference products at inspection. It is also not Java equivalence, source-default execution, clinical validation, or a test of MaCRUISEplus.

The follow-up reuses the converged second-pass TOADS/Dura products and runs the CRUISE core with source-selected ACE and TGDM settings (ACE 25; TGDM initial/inner/central/outer 6/4/7/7, with the selected wrapper mapping central iterations to outer iterations), followed by refinement, surface extraction, and the same quality audit. Its outputs remain under ignored `build/invivo-base-2026-10-10-ras-grid/source-selected-core/`.

The RSS watchdog used by the smoke command samples its direct workflow-runner process, not the full process tree. System-wide `memory_pressure` reported 66–69% free during/after the run; the 7.053 GB runner RSS peak must not be interpreted as total process-tree peak RSS.
