# Synthetic PLUS transfer-mode comparison — 2026-10-12

## Result

The exact same naturally converged iteration-10 registered LesionTOADS
checkpoint was composed twice: once with the recommended physical-affine
transfer and once with the selected historical JIST `UNCHANGED` transfer.
The composer now accepts `--dura-transfer-mode` so this comparison is
restartable without rerunning the expensive lesion loop.

The historical transfer mode completes all downstream stages, but does not
improve this synthetic ribbon. Both modes fail the explicit acceptance audit,
and the `jist-unchanged` result is substantially worse on tissue support and
mesh topology:

| Metric | Physical-affine | JIST unchanged |
| --- | ---: | ---: |
| WM inside inner surface | 99.955% | 88.877% |
| GM between inner and outer surfaces | 96.468% | 64.700% |
| CSF outside outer surface | 99.594% | 99.687% |
| Inner outside central nesting voxels | 1 | 14 |
| Mesh components / Euler characteristic | 4 / −8 | 30 / −884 |
| Surface thickness median / p05 / p95 (mm) | 4.200 / 0.150 / 8.479 | 1.814 / 0.158 / 5.842 |
| Level-set subtraction mean thickness (mm) | 8.623 | 6.612 |

This is a downstream synthetic sensitivity check, not a claim that either
transfer mode reproduces the historical complete workflow. The JIST mode
preserves the cropped 166×206×160 Dura array, then assigns the 172×216×164 T1
reference geometry; the physical-affine mode samples into the T1 reference
grid. The observed change therefore reflects both the selected JIST transfer
semantics and the subsequent Python composition under that geometry contract.
The user's note that historical in-vivo examples may have followed the JIST
path remains plausible, but it is not confirmed by this synthetic run.

The Dura affine and T1 affine are aligned in spacing and axis directions, but
physical mapping from target T1 voxel coordinates to the Dura image is an
integer translation of `(1, -5, 7)` voxels, approximately `(0.94, -4.69,
6.56)` mm or 8.12 mm Euclidean. The unchanged JIST path discards that offset
when it assigns the T1 header. This identifies a concrete source of the two
modes' downstream differences on this synthetic input, without establishing
that the same offset occurred in historical in-vivo runs.

The JIST-mode composition finished all six stages from the saved checkpoint
and reported iteration 10 convergence. Peak sampled RSS was **3.243 GB** under
an 8 GB watchdog. The audit remains synthetic workflow evidence, not Java
equivalence or real-human validation.

## Reproduction

```sh
env PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 \
  --interval-seconds 2 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compose_terminal_synthetic_plus.py \
  --loop-manifest build/plus-atlas-consistent-corrected-loop-iterations8to12-2026-10-12/lesion-registered-loop.json \
  --t1 tests/data/synthetic/plus_atlas_consistent/t1.nii.gz \
  --brain-mask tests/data/synthetic/plus_atlas_consistent/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_atlas_consistent/maseg.nii.gz \
  --output-dir build/plus-atlas-consistent-corrected-jist-unchanged-compose-2026-10-12 \
  --dura-transfer-mode jist-unchanged

.venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/plus-atlas-consistent-corrected-jist-unchanged-compose-2026-10-12/synthetic-plus-terminal-surface-composition.json \
  --report docs/reports/plus-atlas-consistent-corrected-jist-unchanged-ribbon-quality-2026-10-12.json
```

The audited JIST-mode result is in the adjacent
[`ribbon-quality` JSON](plus-atlas-consistent-corrected-jist-unchanged-ribbon-quality-2026-10-12.json).
The selected transfer method comparison is recorded separately in the
[held-fixed JIST transfer report](plus-atlas-consistent-corrected-jist-transfer-java-2026-10-12.md).
