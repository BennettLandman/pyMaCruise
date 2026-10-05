# Config-driven workflow CLI — 2026-10-05

## Implemented behavior

`pymacruise workflow-run --config workflow.json [--resume]` now dispatches the
version-1 workflow config to the existing restartable base or plus runner. It
selects the bundled profile atlas descriptor, preserves the resume flag, and
passes plus-specific `PlusPipelineOptions` to the terminal plus pipeline. If
the output directory contains the preflight-only `run.json`, execution moves
that record to `preflight.json` before the runner writes its execution
manifest. Symlinked manifests and an existing `preflight.json` are rejected.

This closes a previously inert configuration field: `plus_options` was
serialized and validated but the full plus runner always constructed default
options. The configured values now reach the MaCRUISEplus file workflow.

## Focused validation

Command:

```bash
PYTHONPATH=src ./.venv/bin/python -m pytest -q --tb=short \
  tests/test_workflow_preflight_cli.py tests/test_plus_full_workflow.py
```

Result: **12 passed, 0 failed, 0 skipped, 6 expected orientation warnings in
97.25s**. The CLI dispatch tests cover both profiles, bundled descriptor
selection, resume forwarding, preflight-manifest preservation, and configured
plus-option forwarding. The existing plus integration case runs the synthetic
folded-cortex workflow through terminal surfaces and thickness.

## Full regression validation

Implementation commit: `4b86a40f3247e09c3a19949f7e0ddd2485dc0247`.

```bash
PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result: **1,195 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
1,083.83s**. Peak sampled RSS was **15.102 GB**; the watchdog did not
terminate the run.

## Scope limit

This validates config parsing, profile dispatch, option routing, and the
underlying synthetic runners in the test suite. It does not yet run a full
synthetic base or plus workflow by invoking `workflow-run` with a real config
in one end-to-end test. P42a remains Started until that config-launched
integration is checked along with clean-install dependency errors and repeat
resume behavior. No Java/source equivalence or human MRI validation is
claimed.
