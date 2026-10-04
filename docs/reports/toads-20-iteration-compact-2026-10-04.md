# Superseded compact TOADS diagnostic with padding 2

> **Not source-parameter-matched.** The command below omitted `--padding` and
> used the helper default of 2. The selected Java comparison uses padding 10.
> Do not use its convergence result as selected-source evidence. The corrected
> padding-10 continuation is in [the 50-ceiling report](toads-50-iteration-compact-2026-10-04.md).

## Run and result

The registered Python TOADS loop ran on the committed compact `(36,45,34)`
synthetic fixture with the helper's nonmatching padding value of 2. It reached
20 iterations without meeting its `1e-4` threshold; terminal ratio was
`0.000530811725`. This run remains documented to explain the parameter audit,
but it is superseded for source-parameter conclusions.

Reproduce from the repository root:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.1 --report-interval-seconds 15 -- \
  ./.venv/bin/python scripts/trace_synthetic_toads_loop.py \
  --fixture tests/data/synthetic/toads_downsampled_vandy_atlas \
  --atlas data/atlases/cruise-atlas-10obj-toads2008.txt \
  --maximum-iterations 20 --padding 2 \
  --output build/java-reference/python-source-vandy-20iter-low-memory.nii.gz
```

The run took **48.567 seconds** for the registered loop, produced a finite
uint8 RAS classification with labels `[1,5,6,14,15,16,17,18,22,24,25]`, and
had output SHA-256
`cebe386000fc059ef993a6deca72955a1fc947dbce8beb40c4bcef7635287e93`.
Watchdog peak sampled RSS was **1.109 GB**, with no termination. A repeated
run matched the full per-iteration energy-ratio trace and output hash exactly;
the second peak was 1.106 GB. Build outputs and logs are ignored artifacts.

## Bounded-memory trace behavior

The diagnostic previously enabled full iteration-history and per-event
recording unconditionally. A 20-iteration attempt with those options enabled
was terminated by an 8 GB watchdog after sampled RSS reached 8.894 GB. The
trace now disables those diagnostic captures by default and offers explicit
`--retain-iteration-history` and `--capture-events` opt-ins. Per-iteration
energy progress is printed without retaining event records. This changes only
the diagnostic capture policy, not production TOADS code.

The helper default had not matched the selected Java padding. The padding-10
comparison and convergence run supersede this diagnostic for source-aligned
claims. Brain-resolution convergence, the source cause of cross-language
divergence, full-workflow numerical acceptance, and real-MRI validation remain
open.
