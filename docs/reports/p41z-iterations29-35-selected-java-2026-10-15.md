# P41z synthetic registered-loop continuation: iterations 29–35 — 2026-10-15

## Result

Resumed the committed 81×101×77 folded PLUS fixture from the verified
iteration-28 checkpoint. Iteration 29 completed at energy ratio
`0.13199986517429352`; a five-pass segment completed iterations 30–34; and a
separate one-pass segment completed iteration 35. All were nonterminal under
the selected `1e-4` stopping threshold. The observed ratios were:

| Iteration | Energy ratio | Stop |
|---:|---:|---|
| 29 | 0.1319998652 | No |
| 30 | 0.1316104978 | No |
| 31 | 0.1327676475 | No |
| 32 | 0.1324348003 | No |
| 33 | 0.1315448284 | No |
| 34 | 0.1312673688 | No |
| 35 | 0.1322520822 | No |

This extended trace confirms persistent oscillation near 0.132, not movement
toward natural convergence. Each saved checkpoint has **9/9 output hashes
verified** and eight finite RAS NIfTIs. The five-iteration segment peaked at
**1.162 GB** sampled RSS under a 4 GB limit; the final iteration peaked at
**0.879 GB**.

Selected Java comparisons reproduce the 28→29 and 34→35 transitions exactly.
For 34→35, the selected shape samples (5,503,849 values), 267,454 thinning
events, refreshed statistics and membership/energy fields, registration
transform and bounds, 254,377 growth events, and 130 relation rows match; the
Java and Python energy ratio is exactly `0.13225208222866058`. The comparator
peaked at 1.084 GB under the 4 GB limit. These are selected-method
comparisons, not a complete Java application or whole-workflow equivalence.

No new surface-quality audit was run for iteration 35. Earlier ribbon and mesh
diagnostics remain separate from the source-agreement decision and are not
claimed to have improved. The full MaCRUISEplus P41z workflow remains
nonterminal.

## Reproduction

The source-checkout continuations used `pymacruise plus-lesion-registered-loop`
with the prior checkpoint, bundled 12-object atlas descriptor,
`--iterations 99`, `--thinning-limit 1000`, and `--max-passes 1`, except the
29→34 segment used `--max-passes 5`. Every command was run under the repository
RSS watchdog with a 4 GB cap. Selected Java comparisons used
`scripts/compare_p41z_registered_thinning_java.py` with the prior checkpoint
as input and the next checkpoint as `--expected-next-checkpoint`.

Machine outputs remain under ignored `build/p41z-case-normalization/`; no
human imaging data were used.
