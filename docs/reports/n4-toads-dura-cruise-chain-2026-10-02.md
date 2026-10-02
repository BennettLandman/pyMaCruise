# Synthetic N4-to-CRUISE base chain (P18cybx)

## Composition

The compact synthetic CLI sequence now begins with the standalone `n4`
command, runs ANTsPy's actual N4 implementation with the selected default
configuration, then passes that corrected T1 to first-pass registered TOADS.
The N4-derived T1 is also used as the Dura reference. The chain continues
through one Dura outer pass, second-pass TOADS, and the base ACE/GVF/nested
TGDM core. The separate core calls cover both the no-cleanup route and the
optional same-grid MaSeg cleanup route.

The synthetic TOADS runs are limited to two iterations each, and the Dura
stage to one outer pass. This is a compact integration fixture, not a claim of
human-sized convergence.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_loop_cli.py::test_dura_stripped_t1_and_mask_feed_second_toads_cli_pass
```

Result: **1 passed, 0 failed, 0 skipped in 87.76 seconds**.

The test executes N4 rather than mocking it, then successfully consumes its
output in first-pass TOADS and Dura. Both downstream core branches finish and
produce finite float32 inner/central/outer level sets on the original RAS T1
grid. Existing focused N4 tests separately compare CLI output to direct API
output and check the selected configuration.

## Remaining boundaries

Historical SPECTRE registration/transfer remains unresolved. This test does
not establish full TOADS convergence, complete MaCRUISE Java-wrapper
equivalence, visualization surfaces or thickness, or actual-MRI acceptance.
