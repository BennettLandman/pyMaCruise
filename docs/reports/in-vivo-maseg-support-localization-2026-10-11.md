# In-vivo MaSeg support localization — 2026-10-11

## Scope

This is a read-only localization of the human MaSeg support gap using the
already-saved RAS MaSeg, source-faithful refined segmentation, inner/outer
level sets, and supplied brain mask. It does not rerun refinement, alter
algorithm settings, or compare the result with an unavailable historical
whole-pipeline output.

## Result

The saved images are on one 208×309×305, 1 mm RAS grid. Using the existing
support definition (`GM`: inner > 0 and outer ≤ 0; `WM`: inner ≤ 0), the input
MaSeg contains **525,899/642,291 supported GM voxels (81.8786%)** and
**367,582/428,303 supported WM voxels (85.8229%)**. After source-faithful
refinement, recomputing support from labels still carrying selected GM/WM
codes yields **525,899/588,915 GM (89.2996%)** and **367,582/379,771 WM
(96.7904%)**.

Those apparent gains are denominator changes, not improved surface placement:
the number and locations of supported GM and WM voxels are **exactly unchanged**
before and after refinement. The source-faithful label edits remove voxels from
the support denominator. Among the original label cohorts, only 588,915/642,291
GM voxels (91.6897%) and 379,771/428,303 WM voxels (88.6688%) still carry one
of their expected tissue-class codes after refinement; supported voxels within
those cohorts remain at the original 81.8786% and 85.8229%. Therefore the
after-map percentages must not be interpreted as the refinement fixing the
anatomical support gap.

The source-faithful MATLAB `maxind` behavior is visible in the label
transitions. Output codes absent from the input are **1: 39,671 voxels, 2:
11,527, 3: 1,574, and 5: 3**; code 4, which already existed, increases by 146.
The historical row-index assignment can generate values that are not MaSeg
tissue IDs, and in this case removes original class membership from voxels.
Outer clearing also removes unsupported GM labels. The detailed transition
table preserves the before/after source and destination codes.

### Where unsupported voxels lie

Before refinement, unsupported GM comprises 116,392 voxels: 105,169 (90.36%)
are beyond the outer field and 11,223 (9.64%) are on the inner side. Unsupported
WM comprises 60,721 voxels, all outside the inner field. No GM voxels fall in a
crossed-field category (`inner ≤ 0` and `outer > 0`).

| Unsupported cohort before refinement | 0–0.2 | 0.2–0.5 | 0.5–1 | 1–2 | 2–5 | >5 |
|---|---:|---:|---:|---:|---:|---:|
| GM on inner side (11,223) | 6,965 | 3,214 | 916 | 110 | 18 | 0 |
| GM beyond outer (105,169) | 56,036 | 15,486 | 14,694 | 12,187 | 6,380 | 386 |
| WM outside inner (60,721) | 10,763 | 25,302 | 19,697 | 3,913 | 1,018 | 28 |

Bins are magnitudes of the saved signed TGDM fields in their nominal millimeter
scale; this analysis does not independently prove exact Euclidean distance for
every voxel. The RAS octant table in the JSON shows a broad distribution. For
example, the largest GM-beyond-outer octant is RAS with 17,445 voxels (16.6% of
that group); unsupported voxels are not confined to one coarse octant. The
octants split left/right at R=0 and anterior/posterior and superior/inferior at
the supplied mask's median RAS coordinates; these are not atlas regions.

Across the 98 populated GM labels, every label is below 97% support before
refinement and 97 remain below 97% when support is recomputed on the edited
map. For the 17 WM labels, 4 are below 97% before and 5 after; the newly
under-target label 4 reflects an increase of 146 voxels and is part of the
source-faithful row-index concern. Per-label denominators, support, and original
cohort retention are in the CSV.

## Reproduction

From the `pymacruise` repository root:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/localize_in_vivo_maseg_support.py \
  --maseg build/invivo-base-2026-10-10-inputs/ras-grid/orig_target_seg_ras_uint16.nii.gz \
  --refined build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/SEG_refine.nii.gz \
  --inner build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/inner.nii.gz \
  --outer build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/outer.nii.gz \
  --brain-mask build/invivo-base-2026-10-10-inputs/ras-grid/T1_bet_mask_ras.nii.gz \
  --output-dir build/invivo-base-2026-10-10-ras-grid/maseg-support-localization-2026-10-11
```

This saved-product analysis ran in 7.27 seconds. It writes ignored build
artifacts only:

- `build/invivo-base-2026-10-10-ras-grid/maseg-support-localization-2026-10-11/support-localization.json` — input hashes, grid, summary support, spatial/field bins, and caveats.
- `build/invivo-base-2026-10-10-ras-grid/maseg-support-localization-2026-10-11/per-label-support.csv` — before and after per-label support plus original-label-cohort retention.
- `build/invivo-base-2026-10-10-ras-grid/maseg-support-localization-2026-10-11/label-transitions.csv` — counts for every changed source-label → refined-label pair.

The report uses the selected lists from `cortical_separation.py` and repeats
the existing in-vivo audit's sign convention. Its support fraction is not Dice;
the user-agreed Dice criterion is a separate comparison to a reference
segmentation, which is not present here. This analysis also does not establish
clinical validity, source equivalence of the full human pipeline, or exact
distance from every unsupported voxel to a continuous surface.
