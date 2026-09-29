# Full pytest report: P39zn volume smoothing-energy composition

- Date: 2026-09-29
- Implementation commit tested: `856e07b2246642b3dad298e43df52cdaa02eeb7b`
- Command: `.venv/bin/pytest -q`
- Result: **488 passed, 2 warnings in 86.61s**
- Warnings: expected non-RAS topology volume and paint-mask inputs were reoriented to RAS in `tests/test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`.

The new focused selected-Java formula fixture compares every interior voxel
and all three class updates, including all 11 smoothing outputs, as exact
float32 bit patterns. It exercises OBJECT, LESION, and BLACKHOLE smoothing
behavior over fixed synthetic inputs and checks that inputs remain unchanged.
It does not test in-place membership writes, atlas-prior additions, inversion,
normalization, or whole segmentation equivalence.
