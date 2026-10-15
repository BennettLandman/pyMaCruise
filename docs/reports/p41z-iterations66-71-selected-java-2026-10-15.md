# P41z synthetic registered-loop continuation: iterations 66–71 — 2026-10-15

## Result

The folded 81×101×77 synthetic PLUS continuation advanced from iteration 65
to iteration 71. The source-selected stop threshold is `1e-4`; none of the
iterations stopped:

| Iteration | Energy ratio | Stop |
|---:|---:|---|
| 66 | 0.1320327371 | No |
| 67 | 0.1302139759 | No |
| 68 | 0.1302479953 | No |
| 69 | 0.1306818277 | No |
| 70 | 0.1306402832 | No |
| 71 | 0.1316065639 | No |

The trajectory remains near 0.13, over 1,300 times the stop threshold. The
66–70 segment peaked at **1.206 GB** sampled RSS under a 4 GB limit; iteration
71 peaked at **0.797 GB**. Both saved checkpoints verify **9/9 output hashes**
and contain **8/8 finite RAS NIfTIs**.

The selected Java comparison reproduces the 70→71 transition exactly. It
matches 5,503,849 shape samples, 267,443 thinning events, refreshed
statistics and membership/energy fields, registration transform and bounds,
254,324 growth events, and 131 relation rows. Java and Python produce energy
ratio `0.1316065639257431`. The comparator peaked at **0.965 GB** under the
4 GB limit. This is selected-method parity on a matched checkpoint, not a
complete Java application or full-workflow equivalence.

No surface-tail audit was run for iteration 71. The PLUS case remains
nonterminal and this does not imply improved ribbon diagnostics.

## Reproduction

The 66→70 continuation used `pymacruise plus-lesion-registered-loop` from the
iteration-65 checkpoint with the bundled 12-object atlas descriptor,
`--iterations 99`, `--thinning-limit 1000`, `--max-passes 5`, and a 4 GB RSS
watchdog. One additional pass produced iteration 71. Validate 70→71 with
`scripts/compare_p41z_registered_thinning_java.py` using the iteration-70
manifest and iteration-71 checkpoint. Machine outputs remain under ignored
`build/p41z-case-normalization/`.
