# P18cn TOADS uncrop buffer validation — 2026-09-30

## Full-suite evidence

- Command: `.venv/bin/python -m pytest -q`
- Tested implementation commit: `bb9ca9577580e61f66a1c1d55afc7122ce34e6df`
- Result: **643 passed, 2 expected RAS reorientation warnings, 221.63 seconds**.
- The warnings originate in `test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`, which supplies LPS-oriented topology and paint-mask inputs and verifies visible conversion to RAS.

## P18cn — selected `uncropAndBuffer` overloads

The Python `uncrop_toads_buffer` helper reproduces the selected
`MultipleToadPreprocess.uncropAndBuffer` packing for scalar `float[][][]`,
label `byte[][][]`, and channel-first `float[][][][]` crops. It fills the
source-sized output with the supplied background, clips padded support at
source edges, and returns x-fastest voxel buffers; multichannel buffers are
channel-major. The input type is retained, including signed byte values.

The committed asymmetric synthetic case includes negative crop bounds, clipped
x/z support, nonzero background values, two distinguishable channels, and
labels spanning the signed-byte range. The selected Java method-excerpt
harness and Python results match exactly in value and flattened order. Focused
validation: `.venv/bin/python -m pytest -q tests/test_toads_uncrop_buffer.py tests/test_toads_crop.py`
— **7 passed**.

The comparison covers buffer assembly only. It does not validate NIfTI
serialization, MIPAV image import, or end-to-end TOADS output equivalence.
Fixture provenance and the pinned Java source excerpt are documented in
[`tests/data/synthetic/toads_uncrop_buffer/README.md`](../../tests/data/synthetic/toads_uncrop_buffer/README.md).
