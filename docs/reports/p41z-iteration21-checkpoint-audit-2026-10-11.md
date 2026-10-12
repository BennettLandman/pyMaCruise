# P41z iteration-21 checkpoint and selected Java comparison — 2026-10-11

## Runnable synthetic pass

Continued the restartable registered LesionTOADS loop from iteration 20 for one
pass on the committed 71×89×67 RAS synthetic fixture. The command completed in
about 2.5 minutes at **0.746 GB peak RSS** under a **4 GB watchdog**. The
checkpoint is complete; all **9/9** output hashes verify, and all **8/8** NIfTI
outputs are finite and RAS-oriented. The iteration remains nonterminal at
energy ratio **0.1329346597**, up slightly from **0.1327980906** at iteration
20 (threshold `0.0001`).

## Selected Java/source checks

The complete iteration-20→21 selected-method comparison passed under a **4 GB
watchdog** at **0.904 GB peak RSS**. Results match selected Java and the full
Python checkpoint:

- Registered shape sampling: **5,503,849** float values, exact.
- Thinning: **267,570** queue events, **21,422,755-byte** serialization,
  exact.
- Post-thinning statistics and all membership/energy fields: exact float32
  bits; both sides have zero NaNs.
- Registration: transform, rotation, objective, and all 13 shape bounds exact.
- Growth: **254,326** events and **17,558,873-byte** serialization exact;
  segmentation, skeleton, availability, ordering, and class volumes match the
  saved checkpoint.
- Relation refresh: **129** rows and **14** diffusion passes exact; relation
  map and membership table match the saved checkpoint.

These are selected-method source references composed on hash-verified fixed
states, not execution of the complete Java application. The reproduction uses
the shared comparator and its iteration-20 checkpoint:

```sh
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_p41z_registered_thinning_java.py \
  build/p41z-case-normalization/iteration20-complete/lesion-registered-loop.json \
  --output-dir build/p41z-case-normalization/iteration20-next-thinning-java \
  --expected-next-checkpoint build/p41z-case-normalization/iteration21-complete/lesion-registered-loop.json
```

## Ribbon-quality result

The nonterminal surface tail completed at **0.658 GB peak RSS**. Quality still
fails the agreed checks:

- GM-between-surfaces support is **92.905%**, below 97%; WM-inside-inner is
  99.435% and CSF-outside-outer is 99.260%.
- All three meshes are watertight and have one connected component, but Euler
  characteristic is **−18**, not the required 2.
- Surface-thickness median/p05/p95 is **6.501/0.236/12.320 mm**; the p05 and
  p95 thresholds fail. Level-set mean thickness is **4.843 mm** and nesting
  passes.

Compared with iteration 20, mesh connectivity improved from two components to
one, while Euler characteristic worsened from −14 to −18, GM support dropped
from 94.425% to 92.905%, and surface-thickness p95 rose from 11.455 to 12.320
mm. Iteration 21 is source-consistent but neither converged nor anatomically
accepted. Detailed metrics are retained in the ignored
`build/p41z-case-normalization/iteration21-surfaces/ribbon-quality.json`.

This is synthetic validation only. It does not establish real-human PLUS
acceptance or complete-workflow equivalence.
