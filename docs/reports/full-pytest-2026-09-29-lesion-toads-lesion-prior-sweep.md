# Full pytest report: P39zr OBJECT + LESION shape-prior sweep

- Date: 2026-09-29
- Implementation commit tested: `e6eb60e8a962786343f498fb17e9b95f66673e43`
- Command: `.venv/bin/pytest -q`
- Result: **497 passed, 2 warnings in 80.45s**
- Warnings: expected non-RAS topology volume and paint-mask inputs were reoriented to RAS in `tests/test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`.

The selected-Java volume formula harness compares class and auxiliary lesion
membership float32 bits at every interior coordinate after source-order
smoothing, class/lesion atlas-prior additions, direct/prior conversion, shape
gating, and normalization. It includes asymmetric similarities and a
zero-shape lesion voxel. The zero-prior path is checked against P39zq. Distance
modulation, black-hole/OUTMASK classes, registration, and full segmentation
remain outside this evidence.
