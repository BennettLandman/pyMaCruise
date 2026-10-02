# Synthetic MaSeg cleanup through CRUISE core (P18cybw)

## Composition

The compact first-TOADS → Dura → second-TOADS CLI case now exercises both
second-pass core branches: without optional MaSeg cleanup and with
brainstem/cerebellum cleanup enabled. For the cleanup run, the test creates a
same-grid synthetic integer MaSeg map and places each of the eight selected
labels (35, 38–41, 71–73) on a positive GM voxel. `toads-cruise-core` then
zeros those labeled voxels before choosing the GM crop and continues through
(6,18) topology correction, base ACE, generalized GVF, nested TGDM, and
inverse crop.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_loop_cli.py::test_dura_stripped_t1_and_mask_feed_second_toads_cli_pass
```

Result: **1 passed, 0 failed, 0 skipped in 67.50 seconds**.

The CLI provenance confirms `remove_cerebellum=true`, removed counts `[8,8,8]`
for GM, WM-fill, and WM-mask, and `removeBrainstemAndCerebellum` as the first
preparation stage. Final inner/central/outer products are finite float32
NIfTIs on the original T1 RAS grid. This is a synthetic Python composition;
the cleanup primitive's arbitrary-grid/JIST transfer limitation remains as
documented in P18cybr.

## Remaining boundaries

This small case limits each TOADS loop to two iterations and Dura to one outer
pass. It does not include historical SPECTRE registration/transfer, N4
execution inside the composition, human-sized convergence, visualization
surfaces or thickness, real MRI validation, or combined Java-wrapper
comparison. Full MaCRUISE equivalence remains open.
