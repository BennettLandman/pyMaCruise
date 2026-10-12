# P41z registered-loop continuation through iteration 15 — 2026-10-11

## Result

The bounded Python continuation from the saved iteration-11 checkpoint
completed iterations 12–15 without reaching the selected `1e-4` convergence
threshold. Together, iterations 8–15 have ratios in the narrow range
**0.13209–0.13396** and show no downward trend:

| Iteration | Energy ratio | Stop |
| ---: | ---: | --- |
| 8 | 0.13208841 | No |
| 9 | 0.13395977 | No |
| 10 | 0.13252455 | No |
| 11 | 0.13305759 | No |
| 12 | 0.13239078 | No |
| 13 | 0.13359831 | No |
| 14 | 0.13303474 | No |
| 15 | 0.13346674 | No |

Iteration 15 remains nonterminal (`stop=false`, no stop reason). The 12–15
segment took **728.39 seconds** and peaked at **1.073 GB sampled RSS** under a
4 GB watchdog. Each iteration followed the same recorded source order:
thinning, statistics and memberships, shape registration, convergence and
iteration-limit checks, growth, then relation refresh.

## Reproduction

The segment starts from the saved iteration-11 loop manifest and uses its
recorded atlas, source-NaN mode, thinning limit, and convergence threshold.
The bounded continuation call was:

```python
persist_plus_lesion_registered_loop_checkpoint(
    "build/p41z-case-normalization/iterations-0008-to-0011/lesion-registered-loop.json",
    atlas_descriptor,
    "build/p41z-case-normalization/iterations-0012-to-0015",
    n_iter_max=99,
    thinning_limit=1000.0,
    source_zero_variance_nan=True,
    max_passes=4,
    convergence_threshold=1e-4,
    correction_polynomial_degree=None,
    initial_correction_degree=0,
)
```

The run artifacts, including per-iteration progress, are under ignored
`build/p41z-case-normalization/iterations-0012-to-0015/`.

## Interpretation and limits

This demonstrates that the saved P41z state continues to make updates but
does not approach the configured convergence threshold over these eight
iterations. It is a concrete convergence/runtime risk for the selected
checkpoint, not a general failure proof for all inputs or settings. The
iterations 12–15 continuation was not independently replayed through the full
Java wrapper; selected-method Java comparisons elsewhere do not establish
whole-loop equivalence. Source-loop comparison, root-cause diagnosis, and
terminal full-workflow acceptance remain open.
