# Full pytest report: P39zv BLACKHOLE distance sweep

- Date: 2026-09-29
- Implementation commit tested: `c606d9c6c253b114d25581829d59dbec1e202a69`
- Command: `.venv/bin/pytest -q`
- Result: **502 passed, 2 warnings in 79.39s**
- Warnings: expected non-RAS topology volume and paint-mask inputs were reoriented to RAS in `tests/test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`.

The selected-Java synthetic volume harness compares all interior class,
lesion, and black-hole membership bits for the BLACKHOLE distance branch,
with GM enabled and disabled. It covers ventricle precedence, zero-distance
spread, independent brainstem scaling, and the selected ventricle asymmetry
that leaves the black-hole membership unscaled. The no-distance path matches
P39zt. This is not whole LesionTOADS equivalence.
