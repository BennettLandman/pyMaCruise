# Full pytest report: P39zo in-place OBJECT membership sweep

- Date: 2026-09-29
- Implementation commit tested: `5af89244adf65ffab25503e245c82184ee222d8c`
- Command: `.venv/bin/pytest -q`
- Result: **492 passed, 2 warnings in 80.64s**
- Warnings: expected non-RAS topology volume and paint-mask inputs were reoriented to RAS in `tests/test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`.

The OBJECT-only formula harness compares each interior voxel after direct
membership inversion and per-voxel normalization against the selected Java
loop, bit-for-bit. It verifies unchanged boundary planes and a reverse-order
control that differs. This covers fixed per-class data terms, six-face
smoothing, direct relations, and positive shape support only; it does not
cover atlas priors, lesion or other class types, registration, or full
LesionTOADS segmentation.
