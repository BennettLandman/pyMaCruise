# Superseded in-vivo BASE run that retained obliquity — 2026-10-12

> **Superseded for the `--no-regrid-inputs` contract.** This experiment ran
> before the user's clarification that no-regrid must discard oblique
> direction cosines while retaining only closest-RAS axis permutations/flips.
> Its run, metrics, and comparison below document the former preserve-obliquity
> behavior only. The corrected orientation-only human run is pending.

The earlier run directory and outputs remain under ignored `build/` for
reproducibility. The corrected run uses a separate output directory and must
not reuse these saved stages or comparisons.

## Scope and result

This is an end-to-end BASE integration run on the supplied human T1, user-provided brain mask, and MaSeg, with `--no-regrid-inputs`. It tests the native-grid option, which preserves the supplied oblique orthogonal RAS grid. The source audit below shows that the historical V3 route reorients through SPECTRE before TOADS, so this option alone does not reproduce the historical preprocessing. The run completed every workflow stage, passed the stage output checks, and stayed below the 16 GB watchdog. It is **not** a source-equivalence or clinical validation pass.

The T1 grid remained `(184, 256, 256)` with its original oblique affine. All 68 manifest output records exist and their SHA-256 values match; all 51 NIfTI outputs are finite and RAS. Peak sampled RSS was **7.945 GB**. TOADS stopped naturally in both passes: iteration 18 at energy ratio `3.8099166e-05`, then iteration 22 at `1.4343048e-05`.

## Ribbon and source comparison

The quality audit reports nested surfaces except for three voxels where the inner level set falls outside the central level set; central remains inside outer. MaSeg support is **84.05% for GM** between inner and outer and **88.03% for WM** inside inner. The meshes have 74 connected components each, including 72 components under 100 vertices. This is executable human-data integration, but these support and topology measurements show that ribbon quality remains unresolved.

Held-fixed comparison against the downloaded Singularity outputs confirms the source `Raw/T1` and `Raw/Seg` arrays are identical to the supplied inputs. The source level-set NIfTI headers do not preserve the raw T1 oblique affine, so level-set Dice is compared by voxel index, not as a shared-affine physical-space comparison. The index-based nonpositive Dice is **0.97236 inner, 0.97102 central, and 0.97109 outer**; each meets the user’s 0.97 numeric threshold for this comparison. The corresponding symmetric mean surface distances and sampled Hausdorff distances are:

| Surface | Mean nearest-surface distance | Sampled Hausdorff | Result |
|---|---:|---:|---|
| Inner | 0.157238 mm | 4.486190 mm | Fails both position limits |
| Central | 0.241563 mm | 5.832957 mm | Fails both position limits |
| Outer | 0.321702 mm | 7.059284 mm | Fails both position limits |

The held-fixed JIST `IsoSurfaceOnGrid` comparison on the Python level sets passes the surface limits: Java and Python meshes have identical topology and counts; symmetric mean distances are approximately `1.0e-7 mm`, and sampled Hausdorff is `0.00034527 mm` for each surface. This isolates the remaining position mismatch to upstream level-set differences, rather than the Python mesh extraction on these fields. A separate audit found that source TOADS/Dura intermediate NIfTIs use LPS or axis-aligned headers that do not preserve the raw T1 oblique affine; see the [intermediate geometry audit](in-vivo-source-intermediate-geometry-audit-2026-10-12.md). These missing physical mappings limit how strongly voxel-index comparisons can be interpreted.

Refined-label comparison remains below target: 87 of 139 union labels have Dice below 0.97, mean Dice is 0.91958, and labels 1, 2, 3, 5, and 6 appear only in the Python result. The source output has no run logs or manifest, so exact historical executable settings remain unverified. The no-regrid run improves level-set Dice relative to the previously recorded regridded comparison (0.943–0.957), but does not meet the surface-position or refined-segmentation gates.

## Reproduction

Run from the `pymacruise` repository root. The generated human products stay in ignored `build/`; no subject data were copied into the repository.

```sh
PYTHONPATH=src PYTHONNOUSERSITE=1 MPLCONFIGDIR=build/invivo-base-oblique-no-regrid-2026-10-12/mpl \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 2 --report-interval-seconds 45 -- \
  .venv/bin/python -m pymacruise.cli base-run \
  --t1 /Users/landmaba/dev/mipav/sample-macruise/input/T1.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/input/T1_bet_mask.nii.gz \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/input/orig_target_seg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/invivo-base-oblique-no-regrid-2026-10-12 \
  --no-regrid-inputs --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
  --toads-iterations 50 --dura-iterations 5 --ace-diffusion-iterations 25 \
  --gvf-iterations 200 --tgdm-initial-iterations 6 \
  --tgdm-inner-iterations 4 --tgdm-central-iterations 7 \
  --wm-iso-level 0.5 --ace-use-updated-csf-in-ace 2
```

The quality audit and source comparison are saved beside `run.json` under the ignored run directory. The JIST fixed-field comparison is in `jist-surface-comparison/`. The current default remains orthogonal RAS regridding on; this opt-out preserves an oblique orthogonal RAS grid without interpolation. Sheared grids remain unsupported.

The audits can be repeated from the repository root:

```sh
PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/audit_in_vivo_base_products.py \
  build/invivo-base-oblique-no-regrid-2026-10-12/quality-audit-view \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/input/orig_target_seg.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/input/T1_bet_mask.nii.gz \
  --report build/invivo-base-oblique-no-regrid-2026-10-12/quality-audit.json

PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_human_source_output.py \
  --python-output build/invivo-base-oblique-no-regrid-2026-10-12/source-compare-python \
  --python-core build/invivo-base-oblique-no-regrid-2026-10-12/stages/cruise-core/attempt-01 \
  --python-preparation build/invivo-base-oblique-no-regrid-2026-10-12/stages/base-tissue-preparation/attempt-01 \
  --report build/invivo-base-oblique-no-regrid-2026-10-12/source-comparison.json

PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_in_vivo_jist_surfaces.py \
  --core-dir build/invivo-base-oblique-no-regrid-2026-10-12/stages/cruise-core/attempt-01 \
  --surfaces-dir build/invivo-base-oblique-no-regrid-2026-10-12/stages/surface-thickness/attempt-01 \
  --build-dir build/invivo-base-oblique-no-regrid-2026-10-12/jist-surface-comparison \
  --report-path build/invivo-base-oblique-no-regrid-2026-10-12/jist-surface-comparison.md
```
