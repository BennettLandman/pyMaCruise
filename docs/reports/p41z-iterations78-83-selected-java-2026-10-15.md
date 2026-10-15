# P41z synthetic registered-loop continuation: iterations 78–83 — 2026-10-15

## Result

The folded 81×101×77 synthetic PLUS case resumed from iteration 77 and
completed iteration 83. The selected stop threshold remains `1e-4`; none of
iterations 78–83 stopped:

| Iteration | Energy ratio | Stop |
|---:|---:|---|
| 78 | 0.1305855662 | No |
| 79 | 0.1308647245 | No |
| 80 | 0.1308888644 | No |
| 81 | 0.1319803149 | No |
| 82 | 0.1319503486 | No |
| 83 | 0.1302084476 | No |

Iteration 83 remains about 1,302 times above the threshold. The 78–82 segment
peaked at **1.173 GB** sampled RSS under a 4 GB limit; iteration 83 peaked at
**0.794 GB**. Both saved checkpoints verify **9/9 output hashes** and contain
**8/8 finite RAS NIfTIs**.

The selected Java comparison reproduces the 82→83 transition exactly. It
matches 5,503,849 shape samples, 267,570 thinning events, refreshed
statistics and membership/energy, registration transform and bounds, 254,465
growth events, and 130 relation rows. Both sides produce energy ratio
`0.13020844757556915`. The comparator peaked at **0.995 GB** under the 4 GB
limit. This is selected-method parity on a matched checkpoint, not a complete
Java application or whole-workflow equivalence.

No surface-tail audit was run for iteration 83. The complete P41z workflow
remains nonterminal.

## Reproduction

Continue the iteration-77 checkpoint with `pymacruise
plus-lesion-registered-loop`, the bundled 12-object atlas descriptor,
`--iterations 99`, `--thinning-limit 1000`, `--max-passes 5`, and a 4 GB RSS
watchdog to iteration 82. Run one more pass to iteration 83. Validate the
82→83 transition with `scripts/compare_p41z_registered_thinning_java.py` and
the iteration-82/83 manifests. Machine outputs remain under ignored
`build/p41z-case-normalization/`.
