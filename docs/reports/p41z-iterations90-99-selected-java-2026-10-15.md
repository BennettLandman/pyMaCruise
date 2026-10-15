# P41z synthetic registered-loop continuation: iterations 90–99 — 2026-10-15

## Result

The folded 81×101×77 synthetic PLUS continuation reached the configured
99-iteration limit. The final energy ratio was `0.13026626408100128`, far
above the source-selected `1e-4` natural-convergence threshold. The actual
limit-99 checkpoint records `stop_reason=iteration_limit`,
`iteration_limit_stop=true`, and `convergence_stop=false`. This is a clear
nonconvergence result for this selected synthetic state, not a claim that the
algorithm is unportable or that a complete original PLUS run has the same
trajectory.

| Iteration | Energy ratio | Stop |
|---:|---:|---|
| 90 | 0.1310552061 | No |
| 91 | 0.1299555749 | No |
| 92 | 0.1299000680 | No |
| 93 | 0.1303516030 | No |
| 94 | 0.1304168999 | No |
| 95 | 0.1310453415 | No |
| 96 | 0.1311262399 | No |
| 97 | 0.1307514310 | No |
| 98 | 0.1307056844 | No |
| 99 | 0.1302662641 | Iteration cap |

The 90–94 segment peaked at **1.200 GB** sampled RSS under a 4 GB limit; the
95–98 segment peaked at **1.141 GB**; iteration 99 peaked at **0.794 GB**.
All saved checkpoints verify their hashes and contain finite RAS NIfTIs: 9/9
and 8/8 for continuing checkpoints, 7/7 and 6/6 for the terminal product.

## Selected-method comparison at the limit

The terminal output contains no skeleton, availability, ordering, or relation
map, so the existing transition comparator cannot use it as a full
registered-loop checkpoint. To compare the 98→99 update, the same iteration
was rerun from the identical iteration-98 checkpoint with an iteration-100
ceiling, exposing the internal continuation state at iteration 99. This is a
method-level probe, not the actual selected cap-99 run. Both cap-99 and probe
runs produce the exact same iteration-99 energy ratio
`0.13026626408100128`.

The selected Java comparison of the 98→99 probe is exact for 5,503,849 shape
samples, 267,472 thinning events, statistics and membership/energy fields,
registration and bounds, 254,382 growth events, and 130 relation rows. Peak
comparator RSS was **0.958 GB** under 4 GB. This is selected-method parity on
matched checkpoint data; it does not compare the complete Java application
or independently validate Java's iteration-limit wrapper branch.

No iteration-99 surface-tail or ribbon audit was run. P41z remains an
unfinished full-workflow validation row, though its current registered-loop
continuation has now provided an explicit iteration-cap nonconvergence
outcome.

## Reproduction

Continue the iteration-89 checkpoint to iteration 94 with
`pymacruise plus-lesion-registered-loop`, the bundled 12-object atlas
descriptor, `--iterations 99`, `--thinning-limit 1000`, and
`--max-passes 5`, under a 4 GB RSS watchdog. Resume from iteration 94 and use
`--max-passes 4` to iteration 98, then one final pass to 99. For the
selected-method probe only, repeat the final pass from the iteration-98
manifest with `--iterations 100`; compare its internal 99 state using
`scripts/compare_p41z_registered_thinning_java.py`. Artifacts remain under
ignored `build/p41z-case-normalization/`.
