# P41z synthetic registered-loop continuation: iterations 84–89 — 2026-10-15

## Result

The folded 81×101×77 synthetic PLUS loop advanced from iteration 83 to 89.
None of iterations 84–89 reached the selected `1e-4` energy-ratio threshold:

| Iteration | Energy ratio | Stop |
|---:|---:|---|
| 84 | 0.1302186996 | No |
| 85 | 0.1318839937 | No |
| 86 | 0.1318850219 | No |
| 87 | 0.1308952421 | No |
| 88 | 0.1308907419 | No |
| 89 | 0.1310504228 | No |

The energy ratio remains about 1,311 times above the threshold. The 84–88
segment peaked at **1.192 GB** sampled RSS under a 4 GB limit; iteration 89
peaked at **0.794 GB**. Both saved checkpoints verify **9/9 output hashes**
and contain **8/8 finite RAS NIfTIs**.

The selected Java comparison reproduces the 88→89 transition exactly. It
matches 5,503,849 shape samples, 267,420 thinning events, refreshed
statistics and membership/energy, registration transform and bounds, 254,332
growth events, and 130 relation rows. Java and Python produce energy ratio
`0.13105042278766632`. The comparator peaked at **0.965 GB** under the 4 GB
limit. This validates selected methods on matched checkpoints, not a complete
Java application or whole-workflow equivalence.

No surface-tail audit was run for iteration 89. P41z remains nonterminal.

## Reproduction

Continue from the iteration-83 manifest using
`pymacruise plus-lesion-registered-loop`, the bundled 12-object atlas
descriptor, `--iterations 99`, `--thinning-limit 1000`, `--max-passes 5`,
and a 4 GB RSS watchdog through iteration 88. One further pass yields 89.
Validate 88→89 with `scripts/compare_p41z_registered_thinning_java.py` using
the iteration-88 manifest and iteration-89 checkpoint. Machine outputs remain
under ignored `build/p41z-case-normalization/`.
