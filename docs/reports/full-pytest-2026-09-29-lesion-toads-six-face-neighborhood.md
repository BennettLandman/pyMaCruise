# Full pytest report: P39zm six-face membership neighborhoods

- Date: 2026-09-29
- Implementation commit tested: `58a3d757e3a8eebb15d7f44cb72f44df416f27d7`
- Command: `.venv/bin/pytest -q`
- Result: **486 passed, 2 warnings in 80.05s**
- Warnings: expected non-RAS topology volume and paint-mask inputs were reoriented to RAS in `tests/test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`.

The new focused test compares every selected-Java neighborhood stencil in a
coordinate-distinct 5×4×6 float32 synthetic volume, including exact face and
x/y/z visitation order. It also checks that the streaming adapter does not
mutate its input arrays and rejects undersized or spatially mismatched inputs.
The Java harness captures the source loop's reads only; it does not instantiate
or compare a LesionTOADS membership update or full segmentation.
