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

The dispatch/preflight suite passed **11 passed, 0 failed, 0 skipped, 6
expected orientation warnings in 1.35s**. It covers both profiles, bundled
descriptor selection, resume forwarding, preflight-manifest preservation,
unknown config fields, and plus-option forwarding.

The config-launched synthetic integration test passed **1 passed in 73.34s**.
It executes the 25³ folded-cortex plus workflow through terminal surfaces and
thickness, starts from a preflight manifest, verifies saved options and RAS
NIfTI outputs, then runs the same config with `--resume` and verifies a
complete terminal manifest again.

## Full regression validation

Workflow implementation commit: `edb7fba9feb5084f9eddc38e5fc982601d7a1d44`.
Installed-wheel CLI verification commit: `1bc2b98c7665fe6b337fd940d440bf25db268784`.

```bash
PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result on `1bc2b98c7665fe6b337fd940d440bf25db268784`: **1,196 passed, 0 failed,
0 skipped, 11 expected RAS warnings in 1,175.06s**. Peak sampled RSS was
**13.989 GB**; the watchdog did not
terminate the run.

The isolated wheel test additionally resolves packaged atlases and invokes
`python -m pymacruise.cli workflow-run --help` with only the installed target
and an ANTs import marker on the path; it passed **1 test in 23.12s**. The
entire workflow has not yet been executed from the installed wheel.

## Scope limit

This validates config parsing, profile dispatch, option routing, and a full
synthetic plus workflow invoked through `workflow-run`, including repeat
resume. The complete integration test used the source checkout and a bundled-
atlas resolver redirected to the committed miniature synthetic atlas; clean
wheel-install execution remains to be checked. No Java/source equivalence or
human MRI validation is claimed.
