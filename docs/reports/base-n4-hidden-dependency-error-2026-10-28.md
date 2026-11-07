# BASE N4 hidden dependency error — October 28, 2026

Adam reported that `base-run` returned to the shell after regridding, with an
empty `stages/n4/attempt-01/stage.log` and no N4 output. His environment has
not been inspected, so its underlying dependency failure remains unconfirmed.

## Reproduction and correction

Blocking the `ants` import reproduces the empty log and absent terminal error
with the previous source. The process exits with status 2, rather than success.
The N4 CLI uses `argparse.error` for unavailable ANTsPy. That raises
`SystemExit`, bypassing BASE's `except Exception` handler while standard error
is still captured in memory. No failed attempt is recorded.

The correction converts nested CLI exits into reported stage failures,
preserves captured error text in `stage.log`, and marks the stage and parent
run failed. Nonzero returned CLI statuses also fail the stage. N4 scientific
parameters, image geometry handling, atlas selection and downstream algorithms
are unchanged. The existing installed package can be diagnosed with direct
`python -m pymacruise.cli n4`; this bypasses BASE's redirected error handling.

## Validation

The new regression failed before the change because stderr was empty. After
the change it checks both an initial attempt and `--resume`: nonzero exit,
visible dependency error, preserved stage log, failed manifest state, no N4
product and no downstream stages.

```bash
PYTHONPATH="$PWD/src" build/plus-tester-venv/bin/python -m pytest \
  tests/test_cli_n4.py tests/test_base_workflow_restartable.py \
  tests/test_base_workflow_options.py -q
```

Result: **9 passed, 4 expected orientation warnings, 33.09 seconds** using the
clean tester dependency environment with current local source. This includes
successful N4 CLI/API comparison and a complete synthetic BASE workflow with
restart and corruption recovery. The development environment also passed the
same selection (9 passed, 36.31 seconds). This is focused validation on macOS,
not a reproduction of Adam's Linux installation or a fresh human MRI run.

The fix and troubleshooting documentation are local changes; Adam's existing
wheel has not been updated or remotely installed. He does not need another
PyMaCRUISE wheel to expose the underlying error or repair missing Python
dependencies. See [the diagnostic commands](../testing/python-base-checks.md#troubleshoot-a-silent-stop-at-n4).
