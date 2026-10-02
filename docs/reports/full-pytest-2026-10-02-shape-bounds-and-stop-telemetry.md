# Full synthetic suite after TOADS sampler and terminal telemetry changes

## Result

Command: `./.venv/bin/pytest -q`

- **1,101 passed, 0 failed, 0 skipped**
- **11 expected RAS orientation warnings**
- Duration: **1,636.28 seconds (27:16) **
- Tested implementation commit: `f6f900d1156ef44cf4230d460e665fb88385f91b`

This run includes the bounded NumPy batching of the selected TOADS
transformed-shape support scan, exact scalar/batched interpolation parity,
registered-loop terminal energy/cap reporting, the synthetic registered-loop
fixtures, and the integrated N4 → TOADS → Dura → TOADS → CRUISE synthetic
chain. Focused results and the explicit 50-iteration cap outcome are recorded
in [the shape-bound report](toads-shape-bounds-vectorization-2026-10-02.md) and
[the terminal-condition report](toads-loop-terminal-stop-report-2026-10-02.md).

The new 50-iteration result is scoped to the compact 13³ fixture: it reached
the configured iteration limit with energy ratio 0.547289431, above the
1e-4 threshold. That is a convergence-risk signal for this deliberately
limited synthetic setup, not a conclusion about human-sized workflow
convergence. The suite does not include the supplied human scan.

## Source-line estimate

The selected TOADS method-representation estimate remains **3,400–4,800 of
13,651 Java source lines (about 25–35%)**. This change cycle adds zero newly
mapped Java method lines: the scan optimization changes the Python execution
strategy for an already represented operation, and terminal-condition output
connects already implemented control fields. The figure is a coarse
method-level representation estimate, not runtime code coverage; see the
[source accounting note](toads-source-line-accounting-2026-10-01.md).
