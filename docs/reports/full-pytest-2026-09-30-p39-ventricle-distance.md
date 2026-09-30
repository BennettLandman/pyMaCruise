# P39zzzb LesionTOADS ventricle distance validation — 2026-09-30

## Full-suite evidence

- Command: `.venv/bin/python -m pytest -q`
- Tested implementation commit: `d5133c4d266ef4d5f4f86d2b7f8533cfdedf95d4`
- Result: **662 passed, 2 expected RAS reorientation warnings, 160.84 seconds**.
- The warnings come from `test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`, which deliberately provides LPS-oriented topology and paint-mask images and verifies conversion to RAS.
- Focused P39zzzb comparison: `.venv/bin/python -m pytest -q tests/test_lesion_ventricle_distance.py` — **2 passed**.

## P39zzzb — selected ventricle-name seed mask and distance

`lesion_toads_ventricle_distance` reproduces the exact name selection in the
historically misspelled source method
`LesionToadSegmentation.computeVentirclsDistance()`: `Ventricles`,
`LeftVentricle`, `RightVentricle`, and `ThirdVentricle`. It reuses
`lesion_toads_distance_to_mask`, the existing selected-Java-matched byte
distance kernel. The source file SHA-256 is
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.

Three synthetic cases compare the full boolean mask and signed-byte distance
map to a pinned Java method-excerpt harness: all four exact names plus
punctuation/case controls, multiple border/interior seeds, an empty mask, a
typical cap, and the selected cap-1 extra propagation behavior. NIfTI fixtures
have true RAS axis codes and anisotropic spacing; the source distance is in
six-connected voxel steps and ignores physical spacing. Harness and fixture
hashes and regeneration instructions are in
[`tests/data/synthetic/lesion_ventricle_distance/README.md`](../../tests/data/synthetic/lesion_ventricle_distance/README.md).

This validates the atlas-name adapter for supplied labels and the reused
distance kernel. It does not validate LesionTOADS membership modulation or
full segmentation evolution.
