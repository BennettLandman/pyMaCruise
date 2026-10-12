# P41z iteration-20 selected Java handoff — 2026-10-11

## Result

The hash-verified iteration-19 registered checkpoint was replayed through
selected Java/Python methods for registered shape-prior sampling, standard-13
thinning, post-thinning statistics, and the next membership update. The
comparison passed exactly:

- All **5,503,849** registered shape values matched (maximum absolute error
  0.0).
- The serialized thinning result was byte-identical at **21,406,805 bytes**,
  including **97,532** initial queue entries, **267,251** queue events,
  **3,170** topology candidate events, and zero remaining queue entries.
- All seven post-thinning statistics groups matched float32 bits.
- All next-membership maps, energies, denominator, running objective, total
  objective, and energy ratio matched float32 bits. Both implementations had
  zero NaNs.
- The resulting energy ratio is **0.132683739**, above the configured stop
  threshold `0.0001`; the registered loop is still nonterminal.

The run completed under a **4 GB RSS watchdog** with **0.903 GB peak sampled
RSS**. The checkpoint manifest hash is
`b060756ae8d04803c38941bf84a464dfeea59cc00ef831f5c63738ae24c8cab9`. The
machine-readable result in the ignored build directory includes hashes for
the Java harness and selected source files.

## Scope and limitations

This is synthetic, fixed-state selected-method source-comparison evidence. It
extends the iteration-19 statistics comparison through thinning and
membership for the next cycle, but does not compare next-cycle registration,
growth, relation refresh, or a complete live Java loop. It is not real-human
validation, and it does not resolve the existing ribbon-quality failures.

## Reproduction

```sh
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_p41z_registered_thinning_java.py \
  build/p41z-case-normalization/iterations-0016-to-0019/lesion-registered-loop.json \
  --output-dir build/p41z-case-normalization/iteration19-next-thinning-java
```

The comparator verifies checkpoint hashes before reading the registered
state. Generated input/output binaries and `comparison.json` are retained in
the ignored `build/` directory for inspection.
