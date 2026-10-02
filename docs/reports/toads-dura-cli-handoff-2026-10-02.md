# Registered TOADS output through Dura CLI (P18cybp)

## Ported behavior

Validate the operational command-line handoff from the registered TOADS loop
into Dura removal. A synthetic descriptor uses the source-recognized class
names; `toads-loop --dura-output-dir` emits the selected classification and
four 255-scaled tissue maps, then `dura` consumes the exported WM fill, GM, and
WM mask with the original T1. The Dura CLI writes its selected output set on
the native RAS grid.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_loop_cli.py::test_toads_wrapper_outputs_feed_the_dura_cli_and_match_api
```

Result: **1 passed, 0 failed, 0 skipped in 15.90 seconds**.

All seven saved Dura images match the direct `run_dura_removal` API output
exactly after one outer evolution pass; each has the original T1 shape and RAS
axis codes. This validates CLI/data routing over synthetic inputs, not a
combined Java MaCRUISE run or the second TOADS pass. P18cybp remains Started.
