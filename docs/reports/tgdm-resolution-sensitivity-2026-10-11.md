# TGDM resolution sensitivity on a folded-cortex phantom — 2026-10-11

## Question and result

Does the approximately 2 mm voxel spacing in the P41z case, by itself,
explain its severe ribbon-support failure? A controlled 1 mm to 2 mm
downsampling of the already validated 65³ folded-cortex phantom preserves
high label support and mesh topology, so coarse spacing alone does **not**
explain P41z's 0.1123% GM support and 59.51% WM support. The experiment does
show that the source level-set-subtraction thickness measurement is strongly
resolution-sensitive: its mean changes from 4.643 mm at 1 mm to 2.289 mm at
2 mm, failing the 5 mm ±1 mm phantom criterion. The production triangle
thickness summary remains plausible at 2 mm.

## Held-constant setup

The committed [2 mm input fixture](../../tests/data/synthetic/tgdm_resolution_sensitivity_2mm/README.md)
is every second voxel of the 1 mm
[`labeled_folded_cortex`](../../tests/data/synthetic/labeled_folded_cortex/README.md)
case, with the affine scaled by two. This retains the same physical field of
view and labels. Both runs use the same source-selected TGDM parameter set
(initial smoothing 6, initial-inner 4, central 7, central shrink 1, inner 4,
outer 7; WM iso 0.6; zero GVF) and production JIST-derived isosurfaces,
triangle thickness, and level-set thickness. Outputs were written under
ignored `build/`; only synthetic inputs and this report are committed.

Reproduce the 1 mm reference:

```sh
PYTHONPATH=src .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 4 --interval-seconds 0.2 --report-interval-seconds 10 -- \
  .venv/bin/python scripts/run_labeled_folded_cortex.py \
  --fixture "$PWD/tests/data/synthetic/labeled_folded_cortex" \
  --output-dir "$PWD/build/tgdm-resolution-sensitivity-1mm"
```

Repeat with `--fixture "$PWD/tests/data/synthetic/tgdm_resolution_sensitivity_2mm"`
and a separate output directory for the 2 mm result. The 2 mm runner returns
exit code 1 because one explicit acceptance check fails; its output report is
written before the nonzero exit.

## Metrics

| Metric | 1 mm | 2 mm | 2 mm result |
| --- | ---: | ---: | --- |
| WM inside inner | 100.000% | 100.000% | Pass |
| GM between inner and outer | 99.983% | 99.160% | Pass |
| CSF outside outer | 100.000% | 100.000% | Pass |
| Nested voxel regions | 0 violations | 0 violations | Pass |
| Meshes: connected components / Euler / watertight | 1 / 2 / yes | 1 / 2 / yes | Pass |
| Triangle thickness median / p05 / p95 | 4.690 / 4.294 / 5.001 mm | 4.532 / 4.033 / 4.902 mm | Pass |
| Level-set-subtraction thickness mean | 4.643 mm | 2.289 mm | **Fail** |
| Peak sampled RSS / watchdog | 0.379 / 4 GB | 0.187 / 4 GB | Pass |

Exact 1 mm output metrics are in the prior
[labeled folded-cortex report](labeled-folded-cortex-validation-2026-10-09.md).
The complete 2 mm metrics are reproduced in the committed fixture's runner
output when the command above is run; fixture image hashes and geometry are
in `provenance.json`.

## Evidence boundary and next action

This is synthetic Python workflow evidence, not a held-fixed Java comparison
on the downsampled case, whole-workflow equivalence, or real-MRI validation.
The 1 mm case has separate selected-Java gates. P41z's support loss remains
unexplained by this experiment; next, compare the iteration-11 fused seed,
TGDM source settings, and anatomy represented by the coarser crop against a
source-faithful full-resolution case. Do not relax support or thickness gates.
