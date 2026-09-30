# Done-row full-suite evidence audit — 2026-09-30

## Scope and method

This is a static inventory of the roadmap's rows currently marked `Done`. It
does not rerun pytest, judge the technical adequacy of each row's validation,
or change any row status. I inspected the `Done` rows in [`plan.md`](../../plan.md)
and the linked full-suite records under `docs/reports/`.

For a row to count as having unambiguous evidence here, its own plan entry or a
row-specific linked report must identify (1) a complete configured pytest
suite, (2) its numeric pass count, and (3) the exact 40-character implementation
commit tested. Focused test counts do not qualify as the full-suite count.
Abbreviated hashes, suite counts not tied to a particular row/run, and reports
that are only inventories or equivalence analyses do not qualify. This checks
the evidence citation only; it does not verify the historical checkout or
independently establish that a cited command completed successfully.

## Inventory result

| Finding | Count |
|---|---:|
| Rows marked `Done` in `plan.md` | 267 |
| Rows with a full-suite pass count and exact tested commit hash in-row or in an eligible linked report | 267 |
| Rows missing an unambiguous exact-hash citation | 0 |

Most qualifying evidence is recorded inline in the relevant plan entry. For
report-backed entries, the report must identify the plan row and bind its count
and full hash to the same complete-suite run. The principal multi-row evidence
record is [`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md), whose
“Evidence by Done row” table provides row-specific count/hash pairs for P01–P17,
P18a–P18ar, P19–P33, P38, and P39a–P39c. P39a is explicitly identified there as
historical/superseded and is not treated as current Done evidence if its plan
status is reopened. Later row-specific runs are recorded in
[`full-pytest-2026-09-29.md`](full-pytest-2026-09-29.md) (including P18bd) and
[`full-pytest-2026-09-30-lesion-atlas-next-iteration.md`](full-pytest-2026-09-30-lesion-atlas-next-iteration.md)
(P18czl and the other named slices in that run). Other qualifying entries
carry their own complete-suite count and 40-character hash in the plan row.

## Resolution and remaining findings

The initial inventory found P18bh's tested hash abbreviated as `b005798`. The
prefix resolves uniquely in the repository history to
`b00579842e96ae782a1bcc29e2641ced6318e4cf`, an existing commit object. I updated
both the P18bh plan row and its dedicated [full-suite report](full-pytest-2026-09-29-toads-classification-labels.md)
with that exact SHA. The report records **535 passed** in the complete suite.

After this citation repair, all 267 current `Done` rows have an unambiguous
full-suite pass count and exact tested commit hash in-row or in an eligible
row-specific report. This remains a static citation audit; it does not rerun or
independently validate historical pytest runs.

## Evidence records consulted

- [`plan.md`](../../plan.md), all 267 rows currently marked `Done`.
- [`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md), row-specific evidence table.
- [`full-pytest-2026-09-29.md`](full-pytest-2026-09-29.md), separately labeled complete-suite runs.
- [`full-pytest-2026-09-29-toads-classification-labels.md`](full-pytest-2026-09-29-toads-classification-labels.md), P18bh evidence.
- [`full-pytest-2026-09-30-lesion-atlas-next-iteration.md`](full-pytest-2026-09-30-lesion-atlas-next-iteration.md), latest cited lesion-atlas run.

No tests were run for this audit.
