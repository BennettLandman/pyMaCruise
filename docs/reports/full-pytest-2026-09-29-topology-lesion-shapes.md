# Full pytest evidence — 2026-09-29 topology and lesion shape loading

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **385 passed in 83.56 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `ddba29215b37b6fb9b72797a44e83c8d7d0f8f42` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run includes source-derived topology queue, connectivity, threshold-stop,
and float32 corrections with an exact 9×9×9 Java harness comparison. The
comparison report records exact image and queue-event results. It also includes
the Dura WM adapter test and P39e shape-prior NIfTI loader cases for declared
structure order, values, aliases, geometry, RAS reorientation, dimensions, and
invalid/missing paths.

The two warnings are expected: the synthetic Dura adapter case deliberately
supplies LPS-oriented volume and paint images, and checks that both are
reoriented to RAS with a warning. No clinical or user imaging data was used.

This test run establishes test-suite health for the tested commit. It does not
complete P25's broader topology comparison matrix, P18bf's two-call Java
wrapper comparison, or LesionTOADS registration/transformed-prior sampling.
