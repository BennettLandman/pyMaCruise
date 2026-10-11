# P41z iteration-5→6 core replay — 2026-10-11

## Result

The saved synthetic iteration-5 state was replayed through thinning, the
post-thinning statistics refresh, and the GM-enabled membership/objective
update. Each selected Java excerpt matched the Python implementation exactly.
The thinning replay matched the complete serialized initialization and queue
event trace, not only the final labels. Its final segmentation, skeleton,
availability, and ordering also matched the committed state used by the next
two comparisons.

| Stage | Java/Python result |
|---|---|
| Iteration-6 thinning | Exact serialized trace; 98,325 initialization entries and 267,638 queue events |
| Statistics refresh | Exact float32 bits for all seven outputs; 0 NaNs on either side |
| Membership update and objective | Exact float32 bits for every compared output; 0 NaNs on either side |
| Energy ratio | Exact bits: Java and Python `1040630131` (`0.13159732520580292`) |

This confirms the iteration-6 core transition from the supplied checkpoint
under the selected method excerpts. It does not independently validate the
iteration-5 objective: the ratio's previous energy, `3727.62158203125`, comes
from the saved Python iteration-5 summary. Registration, growth, and relation
refresh were not included, and this is not a complete live Java
`AlgorithmLesionToads` run.

## Reproduction

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/replay_p41z_iteration5_to6_core.py
```

Result: success; peak sampled RSS **0.640 GB** under the **4 GB** watchdog,
not terminated for memory. This was a focused source-comparison replay, not a
pytest run or full suite. The implementation baseline was commit
`3465296ede9d62cadaa9863aa8c5ae3676f9ba71`; the standalone replay script was
uncommitted at report time (SHA-256
`1b28ea10706587ce96a0dbc6e1e79cd1c4f9512c5be1fafc207d41b867ea231c`). The
machine-readable evidence is
[`p41z-iteration5-to6-core-replay-2026-10-11.json`](p41z-iteration5-to6-core-replay-2026-10-11.json).

## Input provenance and limits

The replay uses committed RAS NIfTI fixtures under
`tests/data/synthetic/plus_iteration_6_thinning/` and
`tests/data/synthetic/p41z_iteration6_membership_java/`, rather than ignored
build artifacts. The iteration-5 memberships, segmentation, skeleton, and
summary hashes are recorded in the JSON report. The post-thinning fixture is
checked against a fresh thinning result before downstream comparisons run.

The Java references are selected-method excerpts: thinning uses
`Lesion13ClassNextThinningReferenceHarness`, statistics use
`LesionNextIterationStatisticsReferenceHarness`, and membership/objective
uses `LesionNextIterationMembershipReferenceHarness`. The Python and Java
stages are compared on identical saved inputs. Iteration 5 itself, atlas
registration, growth, relation refresh, and full workflow behavior remain
outside this replay.
