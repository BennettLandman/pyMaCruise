# Full pytest report: P39zp OBJECT shape-prior membership sweep

- Date: 2026-09-29
- Implementation commit tested: `7b05f5c1da239cdfd0553e21e617da4b0d65441c`
- Command: `.venv/bin/pytest -q`
- Result: **494 passed, 2 warnings in 79.30s**
- Warnings: expected non-RAS topology volume and paint-mask inputs were reoriented to RAS in `tests/test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`.

The selected-Java formula harness compares every interior OBJECT membership
bit after in-place smoothing, pairwise shape-prior addition, direct/prior
conversion, and voxel normalization. The synthetic fixture includes mixed
relations and a zero-shape class at one interior voxel. A zero-prior/all-direct
control matches the P39zo sweep. Lesion/black-hole classes, auxiliary maps,
registration, and full segmentation remain outside this evidence.
