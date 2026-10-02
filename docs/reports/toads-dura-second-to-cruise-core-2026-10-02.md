# Synthetic first TOADS through CRUISE core CLI chain (P18cybv)

## Composition

The existing compact synthetic CLI case now proceeds through first-pass
registered TOADS, one Dura-removal pass, second-pass registered TOADS, then
`toads-cruise-core`. The final command consumes the second pass's GM,
WM-fill, and WM-mask outputs and runs the selected GM-derived crop, WM
(6,18) topology correction, base ACE, generalized GVF, nested TGDM, and
inverse crop using CLI algorithm defaults. The compact fixture limits each
TOADS run to two iterations and uses one Dura outer iteration. Cropped core
intermediates are intentionally larger than the source grid when the selected
crop edge pads beyond its boundary; final inner/central/outer level sets are
restored to the source T1 grid.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_loop_cli.py::test_dura_stripped_t1_and_mask_feed_second_toads_cli_pass
```

Result: **1 passed, 0 failed, 0 skipped in 50.62 seconds**.

The test confirms that second-pass TOADS products are accepted by the core CLI;
inner, central, and outer level sets are finite float32 RAS NIfTIs on the
original T1 grid; and ACE/GVF crop products share RAS geometry, with a
three-component GVF. This is a synthetic Python CLI integration result, not a
combined Java-wrapper comparison.

## Remaining boundaries

The sequence uses native RAS as its working reference space. It does not
include historical SPECTRE registration/transfer, the optional MaSeg
brainstem/cerebellum cleanup, N4 correction, final visualization surfaces or
thickness, human-sized default TOADS convergence, or actual imaging-data
validation. Full MaCRUISE equivalence remains open.
