# Atlas-aligned synthetic plus ribbon audit — 2026-10-10

## Result

The terminal registered LesionTOADS checkpoint reached iteration 99 but stopped
because it reached the configured iteration limit. Its energy ratio was about
0.1503, so this is not evidence of natural convergence. The checkpoint was
composed through the existing Dura output, native-space transfer, MaSeg mask,
and default source-selected plus surface stages. The composition records the
`PlusPipelineOptions` and `TGDMBaselineParameters`, RAS grid, input/output
hashes, and a run signature. It was composed from persisted stages rather than
one uninterrupted `plus-full-run` invocation.

The atlas-to-MaSeg fixture mapping now follows the descriptor's named regions:
topology classes 5–6 → CSF 5; class 15 (Cerebrum-GM) → cortical GM 100; classes
25–27 (cerebral WM) → cortical WM 44; classes 14, 16–18, 22, and 24 →
additional brain tissue 75; and background → 0. T1 and FLAIR hashes are
unchanged from the earlier loop inputs (`d86d43a0…c0020` and
`ce4bf866…8482b`), so the saved loop remains tied to the same image data.

This labeled fixture still fails ribbon-quality acceptance. It passes RAS/grid,
required-label presence, and voxel-nesting checks. It fails tissue support,
connected genus-zero topology, and both thickness checks:

| Check | Result | Measured value |
|---|---|---:|
| WM inside inner surface | Fail | 52.47% (minimum 97%) |
| GM between inner and outer surfaces | Fail | 0% (minimum 97%) |
| CSF outside outer surface | Pass | 99.92% |
| Voxel nesting | Pass | 0 inner/central or central/outer violations |
| Mesh closure | Pass | all three surfaces watertight |
| Mesh connectivity/topology | Fail | 3 components; Euler characteristic −46 per surface |
| Production triangle thickness | Fail | median 0.470 mm; p05 0.011 mm; p95 0.810 mm |
| Level-set subtraction thickness | Fail | mean 0.821 mm over 246 crossing voxels |

The required synthetic thickness is nominally 5 mm, with median tolerance 1 mm,
p05 ≥2 mm, p95 ≤8 mm, and level-set mean tolerance 1 mm. This result does not
pass those limits. It is synthetic workflow evidence, not Java equivalence or
real-MRI validation.

## Reproduction

The restartable loop checkpoint and composed output are generated under the
ignored `build/` directory. Recompose and audit them with:

```sh
./.venv/bin/python scripts/compose_terminal_synthetic_plus.py \
  --loop-manifest build/plus-vandy-atlas-aligned-v3/loop-chunks/iteration-0100/lesion-registered-loop.json \
  --t1 tests/data/synthetic/plus_vandy_atlas_aligned/t1.nii.gz \
  --brain-mask tests/data/synthetic/plus_vandy_atlas_aligned/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_vandy_atlas_aligned/maseg.nii.gz \
  --output-dir build/plus-vandy-atlas-aligned-v3/terminal-surfaces-v4
./.venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/plus-vandy-atlas-aligned-v3/terminal-surfaces-v4/synthetic-plus-terminal-surface-composition.json \
  --report docs/reports/plus-aligned-folded-ribbon-quality-2026-10-10.json
```

Focused mapping regression: `tests/test_plus_vandy_atlas_aligned_fixture.py` —
**2 passed in 0.31s**. The existing full-profile folded-cortex TGDM case remains
separate evidence and passes its own ribbon checks and held-fixed Java surface
gates; it does not erase this atlas-aligned plus integration failure.

## Next diagnostic

Trace why selected TGDM outputs on this plus handoff collapse the modeled
cortical ribbon to submillimeter thickness and leave cerebral GM outside the
outer surface. Compare intermediate surfaces and selected Java outputs on the
same saved inputs before changing algorithm behavior or relaxing checks.
