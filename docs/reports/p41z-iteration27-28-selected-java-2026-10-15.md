# P41z legacy restart and iterations 27–28 — 2026-10-15

## Legacy checkpoint restart

The saved iteration-26 checkpoint was created before summaries included an
explicit `nonfinite_membership_energy_values` counter. Its membership and
diagnostic NIfTIs are finite (zero NaN/Inf samples), but the current loader
treated the missing counter as a mismatch and refused to resume it.

The loader now measures the diagnostic array directly. A missing legacy count
is accepted only when that array is finite; unreported nonfinite values still
reject, and a present count must match the measured count. The checkpoint
resumed and produced iteration 27. The regression simulates the old summary
schema and exercises hydration on a hash-consistent synthetic checkpoint.

Focused validation:

```sh
PYTHONPATH=src .venv/bin/pytest -q tests/test_plus_lesion_registered_branch.py --tb=short
```

Result: **12 passed, 0 failed, 0 skipped in 11.20 s**. The full suite was not
rerun because this change is confined to checkpoint hydration and its focused
test module covers that contract.

## Restartable synthetic continuation

Both transitions used the source-selected `1e-4` convergence threshold and a
4 GB RSS watchdog, starting from the preceding verified checkpoint:

```sh
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-case-normalization/iteration26-complete/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-case-normalization/iteration27-complete \
  --iterations 99 --thinning-limit 1000 --max-passes 1

MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-case-normalization/iteration27-complete/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-case-normalization/iteration28-complete \
  --iterations 99 --thinning-limit 1000 --max-passes 1
```

Both completed without hitting the watchdog. Peak sampled RSS was **0.820 GB**
for iteration 27 and **0.848 GB** for iteration 28. Each checkpoint's nine
output hashes verify; both contain eight finite RAS NIfTIs. Their energy ratios
were **0.13237354** and **0.13205983**, respectively, well above `1e-4`; neither
iteration stopped. This is a nonterminal synthetic continuation, not a
completed MaCRUISEplus workflow.

## Selected Java comparison

The iteration-27 registered handoff was compared with selected Java methods,
and its Python transition was checked against the independently produced
iteration-28 checkpoint:

```sh
PYTHONPATH=src .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_p41z_registered_thinning_java.py \
  build/p41z-case-normalization/iteration27-complete/lesion-registered-loop.json \
  --output-dir build/p41z-case-normalization/iteration27-next-java \
  --expected-next-checkpoint build/p41z-case-normalization/iteration28-complete/lesion-registered-loop.json
```

At **0.951 GB peak RSS**, all selected comparisons and checkpoint arrays match
exactly: 5,503,849 shape samples, 267,171 thinning events, float32 statistics
and membership fields, registration transform/rotation/bounds, 253,915 growth
events, and the 130-row relation refresh. The Java/Python next-membership
energy ratio is **0.13205983** on both sides. This validates the selected
iteration methods on this fixed synthetic handoff, not a complete Java
application or full BASE/PLUS equivalence.

## Acceptance interpretation and remaining work

The 2026-10-14 owner clarification makes agreement with the original behavior
and stopping rule the acceptance target, including an iteration-cap stop.
Tissue support, mesh connectivity, Euler characteristic, nesting, and absolute
thickness bounds are diagnostic measurements rather than universal port gates.
This two-iteration comparison confirms selected-method parity while the
trajectory remains nonterminal; it does not determine whether the original
complete application follows the same trajectory or reaches the same stop.

Fresh matched full-workflow Java/MATLAB BASE and PLUS references remain open.
This host lacks MATLAB and a container runtime, and the historical downloaded
outputs do not include a matching input/configuration manifest. The source
runtime blocker and human BASE/PLUS comparisons are tracked under P36, P43, and
P44.
