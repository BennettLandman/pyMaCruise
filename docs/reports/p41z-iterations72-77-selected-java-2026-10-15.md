# P41z synthetic registered-loop continuation: iterations 72–77 — 2026-10-15

## Result

The folded 81×101×77 synthetic PLUS case advanced from iteration 71 to 77.
None of these iterations met the source-selected `1e-4` convergence threshold:

| Iteration | Energy ratio | Stop |
|---:|---:|---|
| 72 | 0.1316085011 | No |
| 73 | 0.1309887916 | No |
| 74 | 0.1310203671 | No |
| 75 | 0.1322095990 | No |
| 76 | 0.1322016269 | No |
| 77 | 0.1305921227 | No |

The final ratio remains about 1,306 times the threshold. The 72–76 segment
peaked at **1.197 GB** sampled RSS under a 4 GB limit; iteration 77 peaked at
**0.848 GB**. Both saved checkpoints verify **9/9 output hashes** and contain
**8/8 finite RAS NIfTIs**.

The selected Java comparison reproduces the 76→77 transition exactly. It
matches 5,503,849 shape samples, 267,596 thinning events, membership/energy
and refreshed statistics, registration transform and support bounds, 254,382
growth events, and 131 relation rows. Java and Python produce energy ratio
`0.13059212267398834`. The comparator peaked at **0.950 GB** under the 4 GB
limit. This is selected-method parity on a matched checkpoint, not a complete
Java application or full-workflow equivalence.

No surface-tail audit was run for iteration 77. P41z remains nonterminal; the
earlier tissue and mesh diagnostics have not been refreshed.

## Reproduction

Continue the iteration-71 checkpoint to iteration 76 with
`pymacruise plus-lesion-registered-loop`, the bundled 12-object atlas
descriptor, `--iterations 99`, `--thinning-limit 1000`, `--max-passes 5`,
and a 4 GB RSS watchdog. One additional pass produces iteration 77. Validate
76→77 with `scripts/compare_p41z_registered_thinning_java.py`, using the
iteration-76 manifest and iteration-77 checkpoint. Machine outputs remain
under ignored `build/p41z-case-normalization/`.
