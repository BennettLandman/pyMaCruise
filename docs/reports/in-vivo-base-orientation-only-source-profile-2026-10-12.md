# Human BASE run with orientation-only no-regrid and source-selected settings — 2026-10-12

## Result and scope

The restartable BASE workflow completed on the supplied human T1, matching
brain mask, and MaSeg using the clarified `--no-regrid-inputs` behavior and
selected workflow settings: WM iso **0.5**, ACE updated-CSF mode **2**, 25 ACE
diffusion iterations, and source-faithful segmentation refinement. All **8/8
stages** completed. Both TOADS passes converged naturally at iterations **21**
and **20**, with terminal energy ratios `4.86711651e-05` and
`1.48473373e-05`; neither reached its iteration cap. Peak sampled RSS was
**7.853 GB** under a 16 GB watchdog. All **68/68** manifest output hashes
verified, and all **51/51** NIfTIs were finite and RAS.

This confirms a functioning human BASE workflow with its selected profile. It
does **not** pass whole-pipeline source or ribbon acceptance. The source and
orientation-only T1/MaSeg voxel arrays are exactly identical, but their
affines differ: no-regrid discards the input oblique directions by contract.
Scalar comparisons therefore use voxel indices. For source surface distances,
the Python mesh is mapped from its RAS output through the Python affine back to
voxel indices, then re-expressed under source `Raw/T1` affine. This enables a
same-index comparison in the source geometry for validation; it does not
restore acquisition-world coordinates to the Python output.

## Findings

| Check | Result | Interpretation |
|---|---:|---|
| TOADS pass 1 / 2 | 21 / 20 iterations; ratios `4.8671e-5` / `1.4847e-5` | Both converge naturally |
| Manifest outputs | 68/68 hashes match; 51/51 NIfTIs finite and RAS | Output integrity passes |
| Peak sampled RSS | 7.853 GB | Below the 16 GB limit |
| Source/Python T1 and MaSeg | Exact voxel-array identity; affine differs | Permits index comparison, not same-affine claims |
| Atlas GM / WM / mask Dice vs source | 0.98908 / 1.00000 / 1.00000 | Held-fixed preparation comparison by index |
| Brain-tissue-mask Dice vs source | 0.99354 | Held-fixed preparation comparison by index |
| Inner / central / outer nonpositive level-set Dice | 0.97240 / 0.97119 / 0.97102 | All meet the user's 0.97 Dice threshold by index |
| Refined segmentation | mean label Dice 0.91878; median 0.95992; minimum 0; 86/139 labels below 0.97 | Fails label-wise segmentation tolerance |
| Surface-position mean / sampled Hausdorff, inner | 0.15653 mm / 5.63056 mm | Fails 0.01 mm / 1 mm limits |
| Surface-position mean / sampled Hausdorff, central | 0.24067 mm / 5.83323 mm | Fails 0.01 mm / 1 mm limits |
| Surface-position mean / sampled Hausdorff, outer | 0.32186 mm / 7.05927 mm | Fails 0.01 mm / 1 mm limits |
| MaSeg GM ribbon support; WM-inside-inner support | 84.04%; 88.00% | Ribbon remains incomplete |
| Mesh connected components | 82 per surface; 80–81 under 100 vertices | Severe fragmentation; zero boundary/nonmanifold edges |
| Level-set nesting | 3 inner-outside-central voxels; 0 central-outside-outer | Small inner/central nesting violation |
| Selected JIST surface-extractor replay | Same topology/counts; mean 1.15e-7–1.29e-7 mm; sampled maxima 0.000345–0.000488 mm | Passes extractor limits on fixed level sets |

The selected JIST extractor agreement isolates mesh generation on these exact
fields and passes the requested surface tolerances. The source-levelset
positions still fail, and refinement plus ribbon support remain unresolved.
Source output files have no run log or machine-readable settings manifest, so
the selected Python profile is reproducible but exact historical executable
settings cannot be confirmed from the downloaded result alone.

## Reproduction

Run from the repository root. Input and downloaded source images remain in the
user-provided sample directories and are not committed.

```sh
RUN=build/invivo-base-orientation-only-source-profile-2026-10-12
export PYTHONPATH=src PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 2 --report-interval-seconds 45 -- \
  .venv/bin/python -m pymacruise.cli base-run \
  --t1 /Users/landmaba/dev/mipav/sample-macruise/input/T1.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/input/T1_bet_mask.nii.gz \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/input/orig_target_seg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir "$RUN" --no-regrid-inputs \
  --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
  --toads-iterations 50 --dura-iterations 5 --ace-diffusion-iterations 25 \
  --ace-use-updated-csf-in-ace 2 --wm-iso-level 0.5 \
  --gvf-iterations 200 --tgdm-initial-iterations 6 \
  --tgdm-inner-iterations 4 --tgdm-central-iterations 7
```

The run manifest, per-stage outputs, quality audit, source comparison, and
extractor replay are retained under ignored `build/`. The source comparison
uses `scripts/compare_human_source_output.py --voxel-index-only`; that mode
requires exact orientation-only T1 and MaSeg array identity before comparing.
The held-fixed selected-JIST replay uses `scripts/compare_in_vivo_jist_surfaces.py`
on the run's `stages/cruise-core/attempt-01` and
`stages/surface-thickness/attempt-01` directories.

## Next work

Resolve the remaining source-refinement and surface-position differences,
investigate the fragmented human ribbon, and keep source extractor agreement
separate from upstream segmentation equivalence. Human PLUS validation still
requires a supplied FLAIR input.
