# P41a: synthetic P39 Dura outputs through MaCRUISEplus surfaces

This integration extends the existing crop-ready plus composition upstream to
the P39 Dura output contract. The test creates four Dura channels by calling
`generate_lesion_dura_removal_outputs` with synthetic class memberships, then
passes the selected GM, WM-fill, and WM-mask channels through transfer and
crop preparation. The resulting cropped atlas/TOADS maps feed the selected
topology, boundary, tissue, ACE, GVF, and TGDM stages. TGDM outputs are then
uncropped to the original RAS grid for inner/central/outer surfaces and
thickness products.

Validation command:

```text
./.venv/bin/pytest -q tests/test_plus_p39_to_outputs_integration.py
```

Result: **1 passed in 20.86s**. The test uses committed synthetic 25³ RAS
NIfTI inputs, verifies Dura channel values and stage handoff, finite RAS
volumes, restored full-grid dimensions, RAS+ millimeter meshes, thickness,
and immutability of the supplied images.

The full suite was run with
`./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`.
It passed **1,154 tests, 0 failed, 0 skipped, 11 expected RAS warnings in
834.21s** on tested implementation commit
`677975369a8a11e7776717eed51a44dcee517c69`. Peak sampled RSS was **15.381 GB**;
the watchdog did not terminate the run. See the [full-suite report](full-pytest-2026-10-04-plus-p39-output-to-surfaces.md).

This is connected synthetic-stage evidence, not an end-to-end plus run. The
atlas/TOADS inputs are precomputed, and no FLAIR registration, atlas
registration, LesionTOADS evolution, restart/manifest behavior, Java
whole-stage comparison, or real MRI is exercised. The existing Lewiner marching-cubes versus
JIST `IsoSurfaceOnGrid` substitution risk still applies.
