# Fine-grid labeled plus ribbon audit with source-matched FMM — 2026-10-10

## Result

I recomposed the iteration-6 state from the committed 81×101×77 RAS fixture
through native Dura transfer, MaSeg separation, and all plus surface stages
after switching TGDM to source-ordered indexed-heap fast marching. The
checkpoint remains nonterminal at iteration 6 (energy ratio `0.131597`).
Labels 0/5/44/75/100 are present, nesting passes with no violations, and CSF
outside the outer surface is 99.977%. The remaining ribbon checks still fail:

| Check | Result | Measurement | Target |
|---|---|---:|---:|
| WM inside inner surface | Fail | 59.81% | ≥97% |
| GM between inner and outer surfaces | Fail | 0% | ≥97% |
| Connected genus-zero meshes | Fail | 2 components, Euler −20 per surface | 1 component, Euler 2 |
| Production triangle thickness | Fail | median 0.372 mm; p05 0.150 mm; p95 0.641 mm | median 5±1 mm; p05 ≥2; p95 ≤8 |
| Level-set subtraction thickness | Fail | mean 0.807 mm over 289 crossing voxels | 5±1 mm |

Thus the held-fixed Java TGDM parity correction does not repair the full
synthetic ribbon. The current evidence points to the synthetic tissue/atlas
handoff and resulting surface topology as a separate integration-quality gap.
This is not a failure of the TGDM Java component gate, which passes on the
held-fixed 73×91×63 crop; see the [component comparison](plus-aligned-fine-grid-tgdm-java-2026-10-10.md).
This audit is synthetic only and is not real-MRI validation.

## Reproduction

```sh
./.venv/bin/python scripts/compose_terminal_synthetic_plus.py \
  --loop-manifest build/plus-vandy-atlas-aligned-fine/loop-chunks/iteration-0006/lesion-registered-loop.json \
  --t1 tests/data/synthetic/plus_vandy_atlas_aligned_fine/t1.nii.gz \
  --brain-mask tests/data/synthetic/plus_vandy_atlas_aligned_fine/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_vandy_atlas_aligned_fine/maseg.nii.gz \
  --output-dir build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006-source-fmm \
  --allow-checkpoint

./.venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006-source-fmm/synthetic-plus-checkpoint-surface-composition.json \
  --allow-checkpoint \
  --report docs/reports/plus-aligned-fine-grid-ribbon-quality-source-fmm-2026-10-10.json
```

The full metrics are in the [machine-readable audit](plus-aligned-fine-grid-ribbon-quality-source-fmm-2026-10-10.json).
