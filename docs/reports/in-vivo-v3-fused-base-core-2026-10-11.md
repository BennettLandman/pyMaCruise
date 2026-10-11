# In-vivo V3 fused base-core checkpoint — 2026-10-11

## Scope and result

This checkpoint resumes from the saved second-pass TOADS maps from the existing
human run. It executes V3 MaSeg-derived tissue preparation/fusion, ACE, GVF,
TGDM, source-faithful MaSeg refinement, and the scientific surface/thickness
tail. It does not rerun N4, registration, or TOADS, and it is not a single
restartable workflow-manifest run. The upstream TOADS replay and its source
comparisons are reported separately in the links below.

The MaSeg used here is the run-manifest input
`build/invivo-base-2026-10-10-inputs/ras-grid/orig_target_seg_ras_uint16.nii.gz`,
not the raw oblique source file. It is a nearest-neighbor RAS regrid to the
same 208×309×305 affine as the supplied T1 and brain mask. The saved second-pass
TOADS maps also use that exact grid. Thus the transfer in this human checkpoint
is an identity-grid transfer; this does not validate nonidentity MIPAV
`ReorientVolume` cases.

The source-selected fusion stage matches the selected Java
`CombineMAandTOADS.execute()` output exactly for all three maps over 3,195,483
voxels per map. Fusion support covers 642,290/642,291 MaSeg GM voxels
(99.9998%) and 428,290/428,303 MaSeg WM voxels (99.9970%). This high fusion
coverage does not carry through to the final ribbon. The terminal surface
regions contain 525,899/642,291 input MaSeg GM voxels (81.8786%) and
367,582/428,303 input WM voxels (85.8229%). Source-faithful refinement leaves
these supported counts unchanged; the higher reported fractions against edited
labels (89.2996% GM, 96.7904% WM) come only from smaller denominators. These
containment fractions are not the user-agreed segmentation Dice against a
historical reference, which is unavailable. See the [support localization
report](in-vivo-maseg-support-localization-2026-10-11.md).

All 23 audited NIfTI files are finite and RAS. The three surfaces are
watertight by edge incidence but each has 158 connected components; the largest
component has Euler characteristic −318. The inner surface has 9 voxels outside
the central surface, while the central surface has none outside the outer
surface. The level-set thickness mean is 4.3581 mm (594,033 crossing voxels);
the mesh closest-point thickness mean is 2.6367 mm. These are diagnostic
measurements, not a Java comparison or a pass against an agreed absolute
thickness range. The selected source parameters were used without changing
thresholds to improve support.

This is real-human integration evidence, not whole-workflow equivalence or
clinical validation. Historical full-run reference outputs were unavailable.

## Reproduction

From the `pymacruise` repository root, the corrected preparation is:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python -m pymacruise.cli base-tissue-preparation \
  --maseg build/invivo-base-2026-10-10-inputs/ras-grid/orig_target_seg_ras_uint16.nii.gz \
  --gm build/invivo-base-2026-10-10-ras-grid/stages/toads-second/attempt-01/dura-inputs/stripped_t1_gm.nii.gz \
  --wm-fill build/invivo-base-2026-10-10-ras-grid/stages/toads-second/attempt-01/dura-inputs/stripped_t1_wmfill.nii.gz \
  --wm-mask build/invivo-base-2026-10-10-ras-grid/stages/toads-second/attempt-01/dura-inputs/stripped_t1_wmmask.nii.gz \
  --edge-width 5 --gm-confident-level 0.001 \
  --output-dir build/invivo-base-2026-10-10-ras-grid/source-selected-maseg-toads-preparation-v3-t1-grid-2026-10-11 \
  --provenance build/invivo-base-2026-10-10-ras-grid/source-selected-maseg-toads-preparation-v3-t1-grid-2026-10-11/provenance.json
```

The source fusion comparison is:

```sh
PYTHONPATH=src .venv/bin/python scripts/compare_base_tissue_preparation_saved_java.py \
  --preparation-dir build/invivo-base-2026-10-10-ras-grid/source-selected-maseg-toads-preparation-v3-t1-grid-2026-10-11 \
  --output-dir build/invivo-base-2026-10-10-ras-grid/java-fusion-v3-t1-grid-2026-10-11
```

The source-selected ACE/GVF/TGDM core ran with a 12 GiB watchdog, ACE 25,
GVF lambda 1.0 / 200 iterations / tolerance 1e-6, and TGDM 6/4/7/7 with WM
isovalue 0.6:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 12 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  .venv/bin/python -m pymacruise.cli toads-cruise-core \
  --gm build/invivo-base-2026-10-10-ras-grid/source-selected-maseg-toads-preparation-v3-t1-grid-2026-10-11/fused_gm.nii.gz \
  --wm-fill build/invivo-base-2026-10-10-ras-grid/source-selected-maseg-toads-preparation-v3-t1-grid-2026-10-11/fused_wmfill.nii.gz \
  --wm-mask build/invivo-base-2026-10-10-ras-grid/source-selected-maseg-toads-preparation-v3-t1-grid-2026-10-11/fused_mask.nii.gz \
  --base-preparation-provenance build/invivo-base-2026-10-10-ras-grid/source-selected-maseg-toads-preparation-v3-t1-grid-2026-10-11/provenance.json \
  --edge-width 5 --ace-diffusion-iterations 25 --gvf-lambda 1.0 \
  --gvf-iterations 200 --gvf-solver-tolerance 1e-6 \
  --tgdm-initial-iterations 6 --tgdm-inner-iterations 4 \
  --tgdm-central-iterations 7 \
  --output-dir build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11 \
  --provenance build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/provenance.json
```

Its output directory is
`build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/`;
peak sampled RSS was 2.923 GiB and the watchdog did not terminate it. The
run uses the three `fused_*.nii.gz` maps and the preparation provenance's
MaSeg-derived crop.

The tail commands were:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python -m pymacruise.cli segmentation-refine \
  --segmentation build/invivo-base-2026-10-10-inputs/ras-grid/orig_target_seg_ras_uint16.nii.gz \
  --inner-levelset build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/inner.nii.gz \
  --outer-levelset build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/outer.nii.gz \
  --output build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/SEG_refine.nii.gz \
  --summary build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/SEG_refine-summary.json \
  --source-faithful

PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python -m pymacruise.cli base-surface-thickness \
  --inner-levelset build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/inner.nii.gz \
  --central-levelset build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/central.nii.gz \
  --outer-levelset build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/outer.nii.gz \
  --output-dir build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/surface-thickness
```

The reproducible audit command is:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/audit_in_vivo_base_products.py \
  build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11 \
  --maseg build/invivo-base-2026-10-10-inputs/ras-grid/orig_target_seg_ras_uint16.nii.gz \
  --refined-segmentation build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/SEG_refine.nii.gz \
  --brain-mask build/invivo-base-2026-10-10-inputs/ras-grid/T1_bet_mask_ras.nii.gz \
  --report build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/human-quality-audit-with-refined.json
```

The detailed machine-readable audit and image products stay in ignored
`build/`. Focused tests for the new preparation and restartable workflow wiring
passed **7 tests in 21.85 s**. No full suite was run in this checkpoint.

## Next work

1. The topology-corrected WM-fill and held-fixed TGDM surfaces now match the
   selected Java source exactly/within the agreed surface gates, including the
   source-preserved small components. The 97% MaSeg support target remains
   unmet; localize that support gap without changing source parameters. See the
   [topology correction comparison](in-vivo-fused-core-topology-held-fixed-java-comparison-2026-10-11.md),
   [seed-component audit](in-vivo-current-seed-components-2026-10-11.md), and
   [TGDM comparison](in-vivo-fused-input-tgdm-held-fixed-java-comparison-2026-10-11.md).
2. Add a restartable human manifest for this corrected preparation/core path
   while reusing hashed upstream TOADS products; require resume/hash checks.
3. Continue nonidentity transfer and crop source comparisons before claiming
   general grid equivalence.

See the [TOADS cycle comparison](in-vivo-toads-border-state-2026-10-11.md),
[prior TGDM stage comparison](in-vivo-tgdm-stage-isolation-2026-10-11.md),
[held-fixed fusion comparison](in-vivo-combine-ma-toads-held-fixed-2026-10-11.md),
and [plan](../../plan.md).
