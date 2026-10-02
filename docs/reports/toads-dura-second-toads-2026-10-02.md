# Synthetic TOADS–Dura–TOADS second-pass handoff (P18cybq)

## Ported behavior

Exercise the selected command path corresponding to the two TOADS calls around
Dura removal: run registered TOADS with the supplied N4 T1 and brain mask;
send its exported GM, WM-fill, and WM-mask maps with the original T1 through
the Dura CLI; then run registered TOADS a second time using Dura's stripped T1
and corrected brain mask. Both TOADS passes also write their selected
classification and Dura-input products. All stages use the shared synthetic
RAS grid and a descriptor with tissue names recognized by the Dura wrapper.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_loop_cli.py::test_dura_stripped_t1_and_mask_feed_second_toads_cli_pass
```

Result: **1 passed, 0 failed, 0 skipped in 33.54 seconds**.

The test verifies that Dura output images are accepted as second-pass inputs,
and that the second classification and GM output are native-grid RAS NIfTIs
with uint8 and float32 storage respectively. The adjacent P18cybp test compares
the seven Dura products exactly against direct `run_dura_removal` output after
one outer pass.

This is a synthetic composition of selected Python CLIs, not execution of the
MATLAB/JIST wrapper. It uses native RAS input in place of the historical
SPECTRE-transferred Dura reference; SPECTRE registration and its inverse
transfer remain an explicit working-space deviation. The second TOADS pass
uses the Dura output mask to meet the current Python CLI's required brain-mask
input contract. It does not validate downstream cerebellum removal, GM crop,
topology, ACE, GVF, or TGDM stages, and is not full base-workflow equivalence.
P18cybq remains Started.
