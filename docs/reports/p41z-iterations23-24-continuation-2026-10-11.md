# P41z iterations 23–24 bounded continuation — 2026-10-11

## Runnable synthetic continuation

Continued the restartable synthetic registered-loop checkpoint from iteration
22 for two passes with the source-selected `1e-4` convergence threshold and a
4 GB watchdog:

```sh
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-case-normalization/iteration22-complete/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-case-normalization/iterations23-24 \
  --iterations 99 --thinning-limit 1000 --max-passes 2
```

Both passes completed in **about 6 minutes** at **1.065 GB peak RSS**. Ratios
were **0.132511586** at iteration 23 and **0.132136837** at iteration 24; both
remain far above `0.0001`, so the result is nonterminal. The iteration-24
manifest hash is
`355ebf2e62a9bdaa84b48f035868dfac3df662869004d1503728d8846f12741b`;
all **9/9** output hashes verify, and all **8/8** NIfTIs are finite and RAS on
the 71×89×67 grid. Stage progress records the expected order through
`computeRelations` on both iterations.

## Surface and ribbon audit

The iteration-24 surface tail used **0.671 GB peak RSS**. It remains nonterminal
and fails the agreed ribbon checks:

- WM/GM/CSF support is **99.416% / 94.455% / 99.266%**; GM is below 97%.
- All three meshes are watertight and connected, but Euler characteristic is
  **−14**, not 2.
- Surface thickness median/p05/p95 is **5.717 / 0.239 / 11.442 mm**; p05 and
  p95 fail. Level-set mean is **4.738 mm**, and nesting passes.

These metrics are close to iteration 22: GM support improves by about 0.026
percentage points and Euler changes from −16 to −14, while thickness bounds
remain essentially unchanged. Iterations 23–24 were not separately replayed
with Java: the two-pass output contains only the terminal checkpoint, while
the existing Java comparator validates one continuation at a time. The latest
selected Java comparison remains the exact iteration-21→22 comparison in the
[iteration-22 audit](p41z-iteration22-checkpoint-audit-2026-10-11.md). This is
synthetic integration only, not real-human PLUS validation.
