# Human BASE run with orientation-only no-regrid, default profile — 2026-10-12

## Result and scope

> This run used WM iso 0.6 and ACE mode 0, the BASE defaults. The selected source profile is WM iso 0.5 and ACE mode 2. The geometry, execution, hash, RAS, TOADS-convergence, and fixed-field JIST findings below remain valid; source-levelset, label-refinement, and anatomy metrics are exploratory. The selected-profile rerun completed; see the [source-profile comparison](in-vivo-base-orientation-only-source-profile-2026-10-12.md) for the plan acceptance evidence.

The restartable BASE workflow completed on the supplied human T1, brain mask,
and MaSeg using the clarified `--no-regrid-inputs` behavior. The run manifest
reports all **8/8 stages complete** and natural TOADS convergence on both
passes (21 and 20 iterations). Peak sampled resident memory was **7.589 GB**
under a 16 GB watchdog. All **68/68** manifest output hashes verified; all
**51/51** NIfTIs were finite and RAS. No subject images were copied into the
repository; run products remain under ignored `build/`.

This is a successful end-to-end execution and a failed anatomical/source-
agreement acceptance. It is not clinical validation. The clarified mode
permutes/flips voxel axes to closest RAS, discards oblique direction cosines,
and does not interpolate. The source and orientation-only T1/MaSeg arrays are
voxel-identical, but their affines differ. Source comparisons below therefore
use voxel indices, and surface distances re-express Python mesh voxel
coordinates under the source `Raw/T1` affine. These measurements do not restore
the physical geometry discarded by the no-regrid option.

## Findings

| Check | Result | Interpretation |
|---|---:|---|
| TOADS pass 1 / 2 | 21 / 20 iterations; both met energy threshold | Natural convergence, no iteration cap |
| Manifest outputs | 68 hashes match; 51 NIfTIs finite and RAS | Output-integrity checks pass |
| Peak sampled RSS | 7.589 GB | Below the 16 GB limit |
| Atlas GM / WM / mask Dice vs source | 0.98908 / 1.00000 / 1.00000 | Held-fixed preparation comparison by index |
| Brain-tissue-mask Dice vs source | 0.99354 | Held-fixed preparation comparison by index |
| Inner / central / outer nonpositive level-set Dice | 0.95561 / 0.96466 / 0.96435 | All fail the 0.97 target |
| Refined segmentation | mean label Dice 0.92328; median 0.95757; minimum 0; 91/138 labels below 0.97 | Fails the label-wise 0.97 target |
| Surface-position mean / sampled Hausdorff, inner | 0.27798 mm / 5.63059 mm | Fails 0.01 mm / 1 mm limits |
| Surface-position mean / sampled Hausdorff, central | 0.29184 mm / 5.83023 mm | Fails 0.01 mm / 1 mm limits |
| Surface-position mean / sampled Hausdorff, outer | 0.38183 mm / 7.05927 mm | Fails 0.01 mm / 1 mm limits |
| MaSeg GM ribbon support; WM-inside-inner support | 84.69%; 85.39% | Ribbon support remains poor |
| Mesh connected components | 82 per surface; 80–81 under 100 vertices | Severe fragmentation; no boundary or nonmanifold edges |
| Level-set nesting | 3 inner-outside-central voxels; 0 central-outside-outer | Minor inner/central nesting violation |
| Selected JIST surface extractor replay | Same topology/counts; mean distances 1.15e-7–1.29e-7 mm; sampled maxima 0.000345–0.000488 mm | Passes held-fixed extractor limits; upstream level sets remain discrepant |

The extractor replay isolates mesh extraction on the exact Python level sets:
it agrees with selected JIST `IsoSurfaceOnGrid` under the user's surface
limits. This does not make the produced human anatomy acceptable; the upstream
level sets and refined labels still miss the source targets. The source output
tree contains no run logs or machine-readable settings, so exact historical
executable settings remain unknown. The saved source level-set affines also do
not encode the acquisition's oblique transform.

## Reproduction

Run from the repository root. Human input paths refer to the user-provided
sample data and are intentionally not committed.

```sh
RUN=build/invivo-base-orientation-only-no-regrid-2026-10-12
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
  --gvf-iterations 200 --tgdm-initial-iterations 6 \
  --tgdm-inner-iterations 4 --tgdm-central-iterations 7
```

The retained manifest, per-stage outputs, quality audit, and detailed source
comparison are in the ignored run directory. The source-comparison script now
has an explicit `--voxel-index-only` mode; it refuses to run unless both
orientation-only T1 and MaSeg arrays are proven identical to source inputs.
Reproduce the source comparison and fixed-field extractor replay with:

```sh
PYOUT="$RUN/source-compare-python"
mkdir -p "$PYOUT/surfaces"
ln -sfn "$(pwd)/$RUN/stages/maseg-refinement/attempt-01/SEG_refine.nii.gz" \
  "$PYOUT/SEG_refine.nii.gz"
for name in inner central outer; do
  ln -sfn "$(pwd)/$RUN/stages/surface-thickness/attempt-01/${name}_surface.vtk" \
    "$PYOUT/surfaces/${name}_surface.vtk"
done
ln -sfn "$(pwd)/$RUN/stages/surface-thickness/attempt-01/thickness_levelset.nii.gz" \
  "$PYOUT/surfaces/thickness_levelset.nii.gz"

.venv/bin/python scripts/compare_human_source_output.py \
  --source-output /Users/landmaba/dev/mipav/sample-macruise/output \
  --python-output "$PYOUT" \
  --python-core "$RUN/stages/cruise-core/attempt-01" \
  --python-preparation "$RUN/stages/base-tissue-preparation/attempt-01" \
  --voxel-index-only \
  --python-input-t1 "$RUN/resources/orientation-only-inputs/T1.nii.gz" \
  --python-input-maseg "$RUN/resources/orientation-only-inputs/MaSeg.nii.gz" \
  --report "$RUN/source-comparison.json"

.venv/bin/python scripts/compare_in_vivo_jist_surfaces.py \
  --core-dir "$RUN/stages/cruise-core/attempt-01" \
  --surfaces-dir "$RUN/stages/surface-thickness/attempt-01" \
  --build-dir "$RUN/jist-surface-comparison" \
  --report-path "$RUN/jist-surface-comparison.md"
```

## Next work

Investigate the source-to-Python level-set differences and their effect on
label refinement and ribbon fragmentation. Keep the successful fixed-field
JIST extractor match separate from whole-pipeline source agreement, and
continue to report the no-regrid physical-coordinate limitation explicitly.
