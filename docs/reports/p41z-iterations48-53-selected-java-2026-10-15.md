# P41z synthetic registered-loop continuation: iterations 48–53 — 2026-10-15

## Result

The folded 81×101×77 synthetic PLUS case continued from iteration 47 through
iteration 53. The selected stop threshold is `1e-4`; all six iterations
remained nonterminal:

| Iteration | Energy ratio | Stop |
|---:|---:|---|
| 48 | 0.1319742352 | No |
| 49 | 0.1323952973 | No |
| 50 | 0.1321686059 | No |
| 51 | 0.1310291588 | No |
| 52 | 0.1309793442 | No |
| 53 | 0.1307533830 | No |

The latest ratio is still about 1,308 times the threshold, although it has
slowly decreased since iteration 47. The 49–52 segment peaked at **1.167 GB**
sampled RSS under a 4 GB limit; iteration 53 peaked at **0.794 GB**. The
iteration-48, 49–52, and 53 checkpoints each verify **9/9 output hashes** and
contain **8/8 finite RAS NIfTIs**.

The selected Java comparison reproduces the 52→53 transition exactly. It
matches 5,503,849 shape samples, 267,511 thinning events, all tested
statistics and membership/energy fields, registration transform and bounds,
254,412 growth events, and 129 relation rows. Java and Python report energy
ratio `0.13075338304042816`. The comparator peaked at **0.962 GB** under the
4 GB limit. This is selected-method parity on matched checkpoints, not a
complete Java application or full-workflow equivalence.

No surface-tail audit was run for iteration 53. This continuation is not
evidence that earlier ribbon-quality diagnostics improved; P41z remains
nonterminal.

## Reproduction

Continue the iteration-48 checkpoint with
`pymacruise plus-lesion-registered-loop`, the bundled 12-object atlas
descriptor, `--iterations 99`, `--thinning-limit 1000`,
`--max-passes 4`, and a 4 GB RSS watchdog to iteration 52. Run one more pass
to iteration 53. Validate the 52→53 transition using
`scripts/compare_p41z_registered_thinning_java.py` with the iteration-52
manifest and iteration-53 checkpoint. Machine outputs remain under ignored
`build/p41z-case-normalization/`.
