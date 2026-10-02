# Terminal TOADS to Dura-wrapper outputs (P18cybo)

## Ported behavior

Carry the terminal exact memberships, hard segmentation, and descriptor class
names into the selected `generateDuraRemovalOutputs` behavior. Restore the
classification with the source background value, pack the four CSF/GM/WM-fill/
WM-mask fields with x-fast uncrop semantics, and apply the wrapper's
unconditional 255 scaling. The registered-loop API exposes both raw products
and final wrapper buffers. The `toads-loop --dura-output-dir` option writes the
classification plus four float32 maps on the original native RAS T1 grid.

This exports inputs for Dura removal; it does not execute the Dura algorithm
or the second TOADS pass.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_loop_cli.py
```

Result: **2 passed, 0 failed, 0 skipped in 15.24 seconds**.

The synthetic CLI test checks all five saved images have native T1 dimensions
and RAS geometry, the classification is uint8, and the four Dura input maps
are float32. The underlying product accumulation and 255-scaling adapters
have separate source-reference tests. Full-suite evidence remains pending.
