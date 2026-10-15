# P41z synthetic registered-loop continuation: iterations 42–47 — 2026-10-15

## Result

The folded 81×101×77 synthetic PLUS continuation reached iteration 47. None
of iterations 42–47 met the source-selected `1e-4` stop threshold:

| Iteration | Energy ratio | Stop |
|---:|---:|---|
| 42 | 0.1322709769 | No |
| 43 | 0.1320878565 | No |
| 44 | 0.1318936348 | No |
| 45 | 0.1322690547 | No |
| 46 | 0.1320970207 | No |
| 47 | 0.1321216077 | No |

The energy ratio remains about 1,320 times the selected threshold. The
iterations 43–46 segment peaked at **1.174 GB** sampled RSS under a 4 GB
limit; iteration 47 peaked at **0.812 GB**. The iteration-42, 43–46, and 47
checkpoints each verify **9/9 output hashes** and contain **8/8 finite RAS
NIfTIs**.

Selected Java comparisons reproduce transitions 41→42 and 46→47 exactly. For
46→47, all tested next-membership values and energy, 5,503,849 shape samples,
267,605 thinning events, statistics, registration and bounds, 254,305 growth
events, and 130 relation rows match. Java and Python produce energy ratio
`0.13212160766124725`. The comparator peaked at **0.958 GB** under the 4 GB
limit. This is selected-method parity on matched checkpoints, not a complete
Java application or full-workflow equivalence.

No surface audit was run for iteration 47. Absolute ribbon and mesh metrics
remain diagnostics; this continuation does not establish that they improved.
P41z remains nonterminal.

## Reproduction

The 42→46 segment used `pymacruise plus-lesion-registered-loop` from the
iteration-42 manifest with the bundled 12-object atlas descriptor,
`--iterations 99`, `--thinning-limit 1000`, `--max-passes 4`, and a 4 GB RSS
watchdog. A one-pass continuation produced iteration 47. Selected Java checks
used `scripts/compare_p41z_registered_thinning_java.py` on 41→42 and 46→47.
Machine outputs remain under ignored `build/p41z-case-normalization/`.
