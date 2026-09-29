# Full pytest report: P39zt BLACKHOLE smoothing sweep

- Date: 2026-09-29
- Implementation commit tested: `089507f45ebc7ff659ee1dd3eb5c88fb6b0330fa`
- Command: `.venv/bin/pytest -q`
- Result: **500 passed, 2 warnings in 81.54s**
- Warnings: expected non-RAS topology volume and paint-mask inputs were reoriented to RAS in `tests/test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`.

The selected-Java formula harness compares all interior class, lesion, and
black-hole membership bits after in-place six-face smoothing, separate
auxiliary data terms, shape gating, relation conversion, and source-ordered
normalization. It covers the BLACKHOLE self-neighbor terms and the combined
auxiliary-before-class denominator addition. Atlas priors, distance modulation,
other class types, and full segmentation remain outside this evidence.
