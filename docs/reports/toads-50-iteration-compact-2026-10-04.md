# Compact source-atlas TOADS through the 50-iteration ceiling

## Result

The selected registered Python TOADS loop ran on the committed `(36,45,34)`
source-atlas synthetic fixture with the matched Java comparison settings,
including **padding 10**. Although configured for up to 50 iterations, it met
its `1e-4` energy-ratio threshold and stopped at **iteration 30**, with ratio
**`6.63225946e-5`**.

This is a compact synthetic Python convergence result with source-matched
parameters. It is not source equivalence: Java/Python objective and ratio
values diverge after iteration 2. It does not validate brain-resolution
convergence, whole-workflow Dice/surface tolerances, or real-MRI behavior.

Reproduce from the repository root:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.1 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/trace_synthetic_toads_loop.py \
  --fixture tests/data/synthetic/toads_downsampled_vandy_atlas \
  --atlas data/atlases/cruise-atlas-10obj-toads2008.txt \
  --maximum-iterations 50 --padding 10 \
  --output build/java-reference/python-source-vandy-50iter-padding10.nii.gz
```

The registered loop took **107.128 seconds**. Watchdog peak sampled RSS was
**1.228 GB**, with no termination. The output is finite uint8 RAS NIfTI with
labels `[1,5,6,14,15,16,17,18,22,24,25]` and SHA-256
`f648dba03c4ceb4a453d49ae122ecef83e82905eacc525ca780f922f6e0985d7`.
Build outputs and logs are ignored artifacts.

The first 12 energy ratios match the earlier source-order Java/Python
comparison run exactly, establishing that this is a continuation of the same
Python trace. Java's corresponding ratios are available through iteration 12;
differences become material after iteration 2. Full per-iteration ratios for
the Python run are:

| Iteration | Energy ratio | Iteration | Energy ratio |
|---:|---:|---:|---:|
| 1 | 1.09092224 | 16 | 0.000330082432 |
| 2 | 0.139480755 | 17 | 0.0478523001 |
| 3 | 0.00422830693 | 18 | 0.0457144529 |
| 4 | 0.131389484 | 19 | 0.00251521589 |
| 5 | 0.0919533819 | 20 | 0.00017921369 |
| 6 | 0.041739393 | 21 | 0.0147435721 |
| 7 | 0.00344306417 | 22 | 0.0151027981 |
| 8 | 0.000194382097 | 23 | 0.00021683189 |
| 9 | 0.0143566113 | 24 | 0.000190723717 |
| 10 | 0.0135548962 | 25 | 0.00169247529 |
| 11 | 0.00013966953 | 26 | 0.00164034695 |
| 12 | 0.000304052664 | 27 | 0.000135328824 |
| 13 | 0.00742342137 | 28 | 0.00127355068 |
| 14 | 0.00742778461 | 29 | 0.000822819653 |
| 15 | 0.000229292345 | 30 | 0.0000663225946 |

## Interpretation and next step

The run shows the selected Python loop can meet the source-configured stopping
criterion on this compact fixture when allowed to run, but its energy ratio is
non-monotonic and its trajectory is not yet source-equivalent. The prior
padding-2 20- and 50-iteration reports are diagnostic-only and explicitly
superseded. The production-CLI `(172,216,164)` run with padding 10 reached its
12-iteration cap at ratio `0.0132785104`; brain-resolution convergence is still
unproven. Next, localize the post-iteration-2 difference, then rerun the
brain-resolution TOADS stage before recomputing the full base route.
