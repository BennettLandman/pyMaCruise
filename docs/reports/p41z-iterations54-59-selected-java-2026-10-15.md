# P41z synthetic registered-loop continuation: iterations 54–59 — 2026-10-15

## Result

The folded 81×101×77 synthetic PLUS case resumed from iteration 53 and reached
iteration 59. Its source-selected stopping threshold is `1e-4`; none of these
iterations stopped:

| Iteration | Energy ratio | Stop |
|---:|---:|---|
| 54 | 0.1307082176 | No |
| 55 | 0.1302770227 | No |
| 56 | 0.1302161217 | No |
| 57 | 0.1317305416 | No |
| 58 | 0.1318182200 | No |
| 59 | 0.1304207295 | No |

The ratio remains roughly 1,304 times the threshold. The iterations 54–58
segment peaked at **1.198 GB** sampled RSS under a 4 GB limit; iteration 59
peaked at **0.830 GB**. Both saved checkpoints verify **9/9 output hashes**
and contain **8/8 finite RAS NIfTIs**.

The selected Java comparison reproduces the 58→59 transition exactly. It
matches 5,503,849 shape samples, 267,499 thinning events, membership/energy
and refreshed statistics, registration and support bounds, 254,390 growth
events, and 131 relation rows. The resulting Java/Python energy ratio is
`0.1304207295179367`. The comparator peaked at **0.940 GB** under the 4 GB
limit. This is selected-method parity on matched checkpoints, not a complete
Java application or full-workflow equivalence.

No surface-tail audit was run for iteration 59. The case remains nonterminal;
this does not imply that the earlier ribbon diagnostics improved.

## Reproduction

The 54→58 continuation used `pymacruise plus-lesion-registered-loop` from the
iteration-53 checkpoint with the bundled 12-object atlas descriptor,
`--iterations 99`, `--thinning-limit 1000`, `--max-passes 5`, and a 4 GB RSS
watchdog. One additional pass produced iteration 59. Validate 58→59 with
`scripts/compare_p41z_registered_thinning_java.py`, passing the iteration-58
manifest and iteration-59 checkpoint. Machine outputs remain under ignored
`build/p41z-case-normalization/`.
