# P41z synthetic registered-loop continuation: iterations 60–65 — 2026-10-15

## Result

The 81×101×77 folded synthetic PLUS state resumed from iteration 59 and
reached iteration 65. All six iterations remained above the source-selected
`1e-4` convergence threshold:

| Iteration | Energy ratio | Stop |
|---:|---:|---|
| 60 | 0.1302670687 | No |
| 61 | 0.1314883977 | No |
| 62 | 0.1315573901 | No |
| 63 | 0.1312792450 | No |
| 64 | 0.1312488467 | No |
| 65 | 0.1320632100 | No |

The iteration-65 ratio is 1,321 times the threshold. The 60–64 segment peaked
at **1.177 GB** sampled RSS under a 4 GB limit; iteration 65 peaked at
**0.884 GB**. The iteration-64 and 65 checkpoints each verify **9/9 output
hashes** and contain **8/8 finite RAS NIfTIs**.

The selected Java comparison reproduces the 64→65 transition exactly. It
matches 5,503,849 shape samples, 267,517 thinning events, refreshed
statistics, membership and energy fields, registration transform and bounds,
254,401 growth events, and 132 relation rows. Both produce energy ratio
`0.13206321001052856`. The comparator peaked at **1.084 GB** under a 4 GB
limit. This validates selected methods on a matched checkpoint, not a complete
Java application or full-workflow equivalence.

No surface audit was run for iteration 65. P41z remains nonterminal; this
iteration provides no new evidence about the prior ribbon-quality diagnostics.

## Reproduction

The 60→64 continuation used `pymacruise plus-lesion-registered-loop` from the
iteration-59 checkpoint with the bundled 12-object atlas descriptor,
`--iterations 99`, `--thinning-limit 1000`, `--max-passes 5`, and a 4 GB RSS
watchdog. One additional pass produced iteration 65. Validate 64→65 with
`scripts/compare_p41z_registered_thinning_java.py` using the iteration-64
manifest and iteration-65 checkpoint. Machine outputs remain under ignored
`build/p41z-case-normalization/`.
