# Full pytest report: P39zq OBJECT + LESION smoothing sweep

- Date: 2026-09-29
- Implementation commit tested: `db5f46743761cbe4feebbfb15343524b0787ae5a`
- Command: `.venv/bin/pytest -q`
- Result: **496 passed, 2 warnings in 81.87s**
- Warnings: expected non-RAS topology volume and paint-mask inputs were reoriented to RAS in `tests/test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`.

The selected-Java formula harness compares all class and auxiliary lesion
membership bits at every interior voxel. The fixture covers two OBJECT classes,
one LESION class, direct/prior relation conversion, six-face class and lesion
smoothing, a zero-shape gate, source-order mutation, and per-voxel
normalization. A reverse-order control differs. Atlas-prior additions,
distance maps, black-hole and OUTMASK branches, registration, and whole
segmentation equivalence are not included.
