# Full pytest report: P39zu BLACKHOLE prior sweep

- Date: 2026-09-29
- Implementation commit tested: `9fe0402f83748b3bd134fbfeee918f9ef4665f38`
- Command: `.venv/bin/pytest -q`
- Result: **501 passed, 2 warnings in 80.53s**
- Warnings: expected non-RAS topology volume and paint-mask inputs were reoriented to RAS in `tests/test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`.

The selected-Java synthetic volume harness compares interior class, lesion,
and black-hole membership float32 bits after in-place smoothing, class-pair
shape priors, lesion and black-hole priors, mixed relation conversion, shape
gating, and normalization. It includes the selected BLACKHOLE-class
contribution to ordinary class priors. The zero-prior path matches P39zt.
BLACKHOLE distance modulation and complete segmentation remain open.
