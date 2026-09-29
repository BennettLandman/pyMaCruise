# Full pytest report: P39zs lesion-distance volume sweep

- Date: 2026-09-29
- Implementation commit tested: `23eb94c768ad778e65004044f790c001124f43a7`
- Command: `.venv/bin/pytest -q`
- Result: **498 passed, 2 warnings in 91.20s**
- Warnings: expected non-RAS topology volume and paint-mask inputs were reoriented to RAS in `tests/test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`.

The selected-Java synthetic volume harness compares every interior class and
lesion membership bit with distance modulation applied after membership
conversion and before normalization. It covers ventricle precedence, GM
fallback with the flag on/off, zero-distance spread, independent brainstem
multiplication, direct/prior lesion relations, and the zero-shape gate. The
no-distance path remains covered by P39zr. Black-hole-specific distance
behavior and complete segmentation are not included.
