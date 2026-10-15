# VTK coordinate precision check — 2026-10-15

## Change

The ASCII legacy-VTK writer now serializes RAS coordinates with 17 significant
digits. The previous 10-digit format did not exactly round-trip the Python
binary64 coordinates. `tests/test_surfaces.py::test_vtk_ascii_surface_roundtrip_and_legacy_lps_frame`
now checks exact vertex-array equality after write/read. The focused surface
module passed: **7 passed in 1.18 seconds**.

## Human held-fixed checks

I regenerated the three production meshes from the saved case-100260 BASE
level sets and compared them to the selected Java
`IsoSurfaceOnGrid.solveOriginal(..., 18, 0, false)` source, using the same
float32 fields and the existing 8 GB RSS watchdog. All three surfaces passed
the user's distance limits. The vertex/face counts and topology counts matched
the selected Java output exactly; symmetric mean vertex-to-triangle distances
were `1.08–1.19e-7 mm`, and sampled Hausdorff estimates were
`0.000345–0.000488 mm`. Peak sampled RSS was `0.943 GB`.

I also repeated the selected Java thickness comparison with the newly written
VTK meshes. It still has sparse nearest-hit sensitivity: `1,962 / 240,378`
central nearest-inner-vertex indices differ, with maximum central projected
thickness difference `2.282 mm` and mean absolute difference `0.001397 mm`.
The previous report recorded `1,956` differing indices. Thus VTK decimal
truncation was not the cause, and this change does not resolve the thickness
interaction. Python still matches selected Java bit-for-bit when both receive
the exact same Java mesh coordinates. The user has not set a thickness-field
acceptance threshold; the approved surface-position gates pass. This remains
a precision-sensitive downstream diagnostic, not evidence of whole-workflow
equivalence.

## Reproduction

```bash
PYTHONPATH=src .venv/bin/pytest -q tests/test_surfaces.py --tb=short
```

Human artifacts are retained in ignored `build/vtk-precision-validation-2026-10-15/`.
The held-fixed Java surface and thickness comparison scripts are
`scripts/compare_in_vivo_jist_surfaces.py` and
`scripts/compare_in_vivo_jist_thickness.py`.
