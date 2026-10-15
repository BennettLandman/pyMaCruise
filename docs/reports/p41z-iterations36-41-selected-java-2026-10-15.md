# P41z synthetic registered-loop continuation: iterations 36–41 — 2026-10-15

## Result

The folded 81×101×77 synthetic PLUS case resumed from its hash-verified
iteration-35 state and completed iterations 36–41. No iteration reached the
selected `1e-4` energy-ratio threshold:

| Iteration | Energy ratio | Stop |
|---:|---:|---|
| 36 | 0.1320009232 | No |
| 37 | 0.1316763461 | No |
| 38 | 0.1314593107 | No |
| 39 | 0.1316124201 | No |
| 40 | 0.1314410120 | No |
| 41 | 0.1324421912 | No |

The trajectory continues to oscillate near 0.132, more than 1,300 times the
threshold. The 37–40 segment peaked at **1.182 GB** sampled RSS under a 4 GB
limit; iteration 41 peaked at **0.859 GB**. Both saved checkpoints verify
**9/9 output hashes** and contain **8/8 finite RAS NIfTIs**.

Selected Java comparisons reproduce transitions 35→36 and 40→41 exactly. At
40→41, all tested next-membership fields and energy, 5,503,849 shape samples,
267,561 thinning events, statistics, registration transform and support
bounds, 254,331 growth events, and 128 relation rows match. Java and Python
produce energy ratio `0.1324421912431717`. The comparator peaked at **1.064
GB** under the 4 GB limit. This is selected-method parity on checkpoints, not
a complete Java application or a full-workflow reference.

No surface-tail audit was run for iteration 41. This continuation does not
change the prior ribbon diagnostics or establish natural convergence. The
complete P41z workflow remains open.

## Reproduction

Run `pymacruise plus-lesion-registered-loop` from each preceding checkpoint
with the bundled 12-object atlas descriptor, `--iterations 99`,
`--thinning-limit 1000`, and a 4 GB RSS watchdog. The 36→40 segment used
`--max-passes 4`; the 40→41 segment used `--max-passes 1`. Compare the
35→36 and 40→41 transitions with
`scripts/compare_p41z_registered_thinning_java.py`, passing the prior
checkpoint and the next checkpoint as `--expected-next-checkpoint`.

Artifacts remain under ignored `build/p41z-case-normalization/`; no human
imaging data were used.
