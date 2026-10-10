# Fine-grid labeled synthetic plus run — 2026-10-10

## Scope and result

This is a source-selected synthetic integration attempt on an 81×101×77 RAS
grid, with approximately 2 mm isotropic spacing. The generator resamples the
bundled LesionTOADS atlas and creates matching T1, FLAIR, user brain mask, and
MaSeg labels 0/5/44/75/100. It is synthetic input only. The registered
LesionTOADS loop reached iteration 6 in restartable chunks and did not
converge (energy ratio **0.131597**, selected threshold **0.0001**). The
checkpointed raw-input workflow ran the source-selected Dura transfer, MaSeg
separation, ACE/GVF/TGDM surface evolution, mesh extraction, and thickness
stages. This did not reach a terminal 99-iteration workflow state.

The explicit ribbon audit did **not** pass:

| Check | Result | Measurement |
| --- | --- | ---: |
| Required MaSeg labels | Pass | 0, 5, 44, 75, 100 present |
| Voxel nesting | Pass | 0 violations |
| CSF outside outer surface | Pass | 99.977% (minimum 97%) |
| WM inside inner surface | Fail | 59.81% (minimum 97%) |
| GM between inner and outer surfaces | Fail | 0% (minimum 97%) |
| Mesh closure | Pass | all three watertight; no boundary or nonmanifold edges |
| Connected genus-zero mesh | Fail | 2 components and Euler characteristic −20 per surface |
| Production triangle thickness | Fail | median 0.372 mm; p05 0.150 mm; p95 0.641 mm |
| Level-set subtraction thickness | Fail | mean 0.807 mm across 289 crossing voxels |

The thickness target remains 5 mm (median within 1 mm; p05 ≥2 mm; p95 ≤8 mm;
level-set mean within 1 mm). Comparing the earlier iteration-3 checkpoint
audit with this iteration-6 audit shows little change: GM support remains 0%,
WM support is 59.7–59.8%, and median thickness is 0.372–0.375 mm. This case
therefore supplies a concrete synthetic integration failure signal. It is not
Java-equivalence evidence or human-MRI validation. The already passing 65³
analytic folded-cortex TGDM fixture remains separate evidence for the isolated
surface stages.

## Restartable run and audit

The first attempt at this finer grid exposed an overstrict guard in source
zero-variance mode: it rejected a later iteration after all current standard
deviations became nonzero. The guard now requires GMflag but does not require a
zero deviation to remain present. The synthetic fixture generator accepts an
output directory and target grid, so the committed case can be rebuilt with:

```sh
./.venv/bin/python tests/data/synthetic/plus_vandy_atlas_aligned/generate_fixture.py \
  --output-dir tests/data/synthetic/plus_vandy_atlas_aligned_fine \
  --shape 81 101 77
```

Run the raw-input plus path to a two-iteration checkpoint:

```sh
./.venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_vandy_atlas_aligned_fine/t1.nii.gz \
  --flair tests/data/synthetic/plus_vandy_atlas_aligned_fine/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_vandy_atlas_aligned_fine/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_vandy_atlas_aligned_fine/maseg.nii.gz \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/plus-vandy-atlas-aligned-fine-2iter \
  --registration-iterations 2 --shrink-factors 1 --smoothing-sigmas 0 \
  --registration-sampling-rate 1 --n4-shrink-factor 1 --n4-iterations 2 \
  --n4-spline-mesh 1 1 1 --crop-padding 1 --init-align-iterations 0 \
  --pyramid-levels 2 --lesion-iterations 1 --first-thinning-limit 1000 \
  --loop-thinning-limit 1000 --max-passes 2
```

The registered-loop CLI extends that saved state in additional bounded chunks.
The iteration-4/5 checkpoint was generated with:

```sh
./.venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/plus-vandy-atlas-aligned-fine-2iter/lesion_workflow/registered-loop/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/plus-vandy-atlas-aligned-fine/loop-chunks/iteration-0005 \
  --iterations 99 --thinning-limit 1000 --max-passes 2
```

Iteration 6 uses that manifest as its input and writes to
`build/plus-vandy-atlas-aligned-fine/loop-chunks/iteration-0006/`:

```sh
./.venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/plus-vandy-atlas-aligned-fine/loop-chunks/iteration-0005/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/plus-vandy-atlas-aligned-fine/loop-chunks/iteration-0006 \
  --iterations 99 --thinning-limit 1000 --max-passes 1
```

Both checkpoints are hash-checked and resumable. The iteration-6 state was
composed through Dura, native transfer, MaSeg, and the plus surface workflow
using:

```sh
./.venv/bin/python scripts/compose_terminal_synthetic_plus.py \
  --loop-manifest build/plus-vandy-atlas-aligned-fine/loop-chunks/iteration-0006/lesion-registered-loop.json \
  --t1 tests/data/synthetic/plus_vandy_atlas_aligned_fine/t1.nii.gz \
  --brain-mask tests/data/synthetic/plus_vandy_atlas_aligned_fine/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_vandy_atlas_aligned_fine/maseg.nii.gz \
  --output-dir build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006 \
  --allow-checkpoint
```

The audit requires the explicit checkpoint opt-in and records nonterminal
status, completed iteration, and loop energy ratio:

```sh
./.venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006/synthetic-plus-checkpoint-surface-composition.json \
  --allow-checkpoint \
  --report docs/reports/plus-aligned-fine-grid-ribbon-quality-2026-10-10.json
```

The run stayed below a **1.1 GB** sampled RSS peak under a 12 GB watchdog. The
separate audit stayed below **0.2 GB** under an 8 GB watchdog. A bounded focused
regression command on implementation commit `072ba03` passed **6 tests**:

```sh
./.venv/bin/python -m pytest tests/test_lesion_standard13_voxel.py -q
```

The full suite was not run. The machine-readable quality audit is
[`plus-aligned-fine-grid-ribbon-quality-2026-10-10.json`](plus-aligned-fine-grid-ribbon-quality-2026-10-10.json).
