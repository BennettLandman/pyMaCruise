# P39zzzc LesionTOADS GM distance validation — 2026-09-30

## Full-suite evidence

- Command: `.venv/bin/python -m pytest -q`
- Tested implementation commit: `4665da22d2c2e30743fb541cb08ae72a01540c0e`
- Result: **664 passed, 2 expected RAS reorientation warnings, 164.97 seconds**.
- The warnings come from `test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`, which deliberately supplies LPS-oriented topology and paint-mask images and verifies conversion to RAS.
- Focused P39zzzc comparison: `.venv/bin/python -m pytest -q tests/test_lesion_gm_distance.py` — **2 passed**.

## P39zzzc — selected GM-name seed mask and distance

`lesion_toads_gm_distance` reproduces the exact atlas-name selection in
`LesionToadSegmentation.computeGMDistance()` for the 13 source names and
reuses the existing selected-Java-matched byte `fastDistance` kernel. The
selected source file SHA-256 is
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.

The committed synthetic RAS NIfTI cases compare every boolean mask and
signed-byte distance voxel with a pinned Java method-excerpt reference.
They cover all exact names, case/punctuation near misses, multiple and
boundary seeds, empty support, cap 1's extra propagation, and a typical cap.
The distance is measured in six-connected voxel steps; as in the selected
Java helper, image spacing is not used. Harness and fixture hashes and
regeneration details are recorded in
[`tests/data/synthetic/lesion_gm_distance/README.md`](../../tests/data/synthetic/lesion_gm_distance/README.md).

This validates the supplied-label GM mask and distance adapter only. It does
not validate the subsequent GM membership attenuation or full LesionTOADS
evolution.
