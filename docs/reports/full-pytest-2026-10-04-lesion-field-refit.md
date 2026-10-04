# Full pytest evaluation after LesionTOADS adaptive field-refit composition

## Result

The repository-wide suite passed against implementation commit
40589ebbf5a567b93d9a7dd44cc86a5f9eb6ec0e:

    1143 passed, 11 warnings in 820.91s (0:13:40)

The warnings are expected RAS reorientation notices from orientation-contract
tests. The suite used only synthetic data and did not run the supplied human
MRI.

Command:

    ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q

The watchdog completed normally with exit code 0. Peak sampled RSS was
14.970 GB; the 24 GB watchdog did not terminate the process. The host maximum
is 40 GB.

## Scope of evidence

This is repository-wide regression evidence after adding descriptor-aware
13-class LesionTOADS exact field fitting and composing its adaptive degree
gate with ordered per-channel updates. The source-formula fixture verifies the
field kernel and branch wiring on synthetic data. The LesionTOADS iterative
workflow, MaCRUISEplus integration, source-equivalent full workflow, and
real-MRI validation remain open.

No newly counted selected base-TOADS Java source lines were added in this
cycle. The reported coarse estimate therefore remains 3,464–4,864 of 13,651
(about 25–36%); its denominator explicitly excludes the separate LesionTOADS
implementation. This estimate is method-level planning, not executed Python
coverage or workflow completion.
