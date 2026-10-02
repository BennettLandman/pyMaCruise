# Standalone registered TOADS loop CLI (P18cyn)

## Ported behavior

Add `pymacruise toads-loop` for the selected single-channel, default-atlas
registered TOADS path. The command accepts an N4-corrected T1, a matching
user-supplied brain mask, an optional descriptor override, and selected
registration/loop limits. It uses the bundled base TOADS atlas by default,
performs atlas registration and the registered loop, restores terminal labels
to the native RAS grid, and writes one compressed NIfTI classification.

This command does not run N4, dura removal, the second MaCRUISE TOADS pass,
or emit the CRUISE products. It is a standalone TOADS interface, not the full
base/plus workflow CLI or workflow-equivalence claim.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_loop_cli.py
```

Result: **2 passed, 0 failed, 0 skipped in 17.81 seconds**.

The synthetic run uses a compact 11-class atlas and one registration level;
the written output has the input T1 shape, `uint8` labels, and RAS orientation.
A separate help check verifies that corrected T1, brain mask, and output are
required arguments. Full-suite and packaged-install CLI evidence is pending.
