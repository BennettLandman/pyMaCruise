# Human BASE MaCRUISE run, RAS-grid attempt 02

## Result

The restartable BASE workflow completed on the supplied human T1, brain mask, and MaSeg. `run.json` records status `complete`; all seven stages completed and passed their stage checks. Both TOADS passes stopped by the energy threshold, not the iteration cap. Peak sampled RSS was **7.302 GB** under a 16 GB watchdog; it did not terminate the run.

This establishes that the current Python BASE route can complete on this human input. It does not establish whole-pipeline Java equivalence, the agreed segmentation Dice or surface-position tolerances, clinical validity, or completion of the plus workflow.

## Reproduction

From `pymacruise/`:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 1 --report-interval-seconds 60 -- \
  .venv/bin/python -m pymacruise.cli base-run \
  --t1 build/invivo-base-source-v3-human-inputs-2026-10-11/ras-grid/T1_raw_ras.nii.gz \
  --brain-mask build/invivo-base-source-v3-human-inputs-2026-10-11/ras-grid/T1_bet_mask_ras.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --maseg build/invivo-base-source-v3-human-inputs-2026-10-11/ras-grid/orig_target_seg_ras_uint16.nii.gz \
  --output-dir build/invivo-base-source-v3-human-rasgrid-2026-10-11 \
  --resume --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
  --toads-iterations 50 --dura-iterations 5 --dura-membership-level 0.3 \
  --crop-edge-width 5 --ace-diffusion-iterations 25 --gvf-iterations 200 \
  --tgdm-initial-iterations 6 --tgdm-inner-iterations 4 \
  --tgdm-central-iterations 7 --wm-iso-level 0.5 \
  --ace-use-updated-csf-in-ace 2 --source-faithful-refinement
```

The manifest records the input hashes, stage attempts, outputs, and RAS convention. The four inputs are the 1-mm RAS T1, matching user-supplied brain mask, RAS MaSeg, and bundled atlas descriptor. T1, mask, and MaSeg are each `208×309×305`; the manifest SHA-256 values are recorded in `run.json`.

## Workflow and source-comparison evidence

- All stages completed on attempt 02: N4, TOADS first pass, Dura, TOADS second pass, base tissue preparation, CRUISE core (ACE/GVF/TGDM), source-faithful MaSeg refinement, and surface thickness.
- First TOADS pass: 24 iterations; terminal ratio `5.83754409e-05`; `energy_threshold_reached=true`, `iteration_limit_reached=false`.
- Second TOADS pass: 26 iterations; terminal ratio `5.61925735e-05`; `energy_threshold_reached=true`, `iteration_limit_reached=false`.
- Selected JIST `ReorientVolume.solve(UNCHANGED, NEAREST_NEIGHBOR)` versus Python transfer on the saved human grids: GM, WM-fill, and WM-mask each compare **19,602,960 voxels**, bitwise equal, zero differing voxels, and zero maximum/mean absolute error. The grids are already aligned RAS grids; this does not test changed orientation, origin, spacing, or shape. Evidence: `build/java-reference/reorient-volume-human-transfer/human-attempt-02/comparison.json`.
- Selected JIST `CubicVolumeCropper` versus Python base-preparation crop: crop bounds match exactly at `[[30,170],[58,230],[105,235]]`, producing `141×173×131`; all seven arrays (three atlas maps, three transferred TOADS maps, and MaSeg) have zero differing voxels. The crop is interior to the image, so boundary padding is not exercised. The harness compares saved Python transfer/separation products; it does not replay those stages. Evidence: `build/java-reference/human-base-cropper/invivo-base-source-v3-human-rasgrid-2026-10-11-attempt-02/comparison.json`.

## Human ribbon-quality audit

The audit checked finite NIfTIs, geometry, nesting, label support, mesh topology, and thickness. All audited images were finite and RAS on the `208×309×305` grid; `bad_nifti=0`. These diagnostics are not Dice scores.

- Input MaSeg support between inner and outer surfaces: **529,145 / 642,291 GM voxels = 82.384%**. Input WM labels inside the inner surface: **378,800 / 428,303 = 88.442%**.
- On the source-refined MaSeg, those support fractions are **89.044% GM** and **96.473% WM**. This is a separate label-support diagnostic, not comparison to ground truth.
- Nesting: zero inner-outside-central and zero central-outside-outer voxels. The inner surface has 305 voxels outside the supplied brain mask.
- Each mesh has zero boundary and nonmanifold edges, but each has **136 connected components** (Euler characteristic `-28`); 135 inner components contain 844 vertices in total. The component count and small fragments remain a ribbon-quality concern.
- Mesh-derived symmetric closest-surface thickness mean: **2.631 mm**. Level-set thickness mean: **4.191 mm**. These are thickness measurements, not mean surface-position error against a reference.
- Dura brain mask Dice with the supplied mask: **0.95198** (37,736 Dura-only and 85,349 supplied-only voxels). This is a mask comparison, not the required cortical segmentation Dice.
- Visual QA overlays show the three surfaces following folded cortical anatomy in sagittal, coronal, and axial views. The separate 3D panels use display-only vertex clustering and do not alter saved surfaces. Files are in `build/qa-human-ribbon-current-attempt-02-2026-10-11/`.

No complete Java/MATLAB output for this same input is available for whole-workflow comparison. The audit therefore records `source_reference_comparison=unavailable`. The agreed minimum Dice `0.97`, Hausdorff `1 mm`, and mean surface distance `0.01 mm` cannot be declared passed from these diagnostics; matching source surfaces/segmentations are still needed.

## Next work

Use the successful restartable case to compare remaining available JIST stages and diagnose the 136-component ribbon. Then run and validate the plus workflow. Keep synthetic, source-equivalence, and human-data evidence separate. The overall port is not complete.
