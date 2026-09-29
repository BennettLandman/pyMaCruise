# Full pytest evidence — 2026-09-29 follow-up

## Combined Dura and LesionTOADS loader slices

| Field | Value |
|---|---|
| Command | `.venv/bin/python -m pytest -q` from the `pymacruise/` repository root |
| Result | **371 passed in 75.25 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `ca391935960baf39c99e880fd83cb5b20b496ebf` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

The run includes the P18be Dura WM topology adapter and the P39d LesionTOADS
topology-atlas NIfTI loader. The Dura test verifies explicit wrapper settings,
RAS warning/reorientation, output dtype, and geometry while invoking the
existing `topology_correct` implementation. It does not compare the corrected
voxels against a Java topology harness. The P39d tests validate bundled
`.nii.gz` loading and RAS behavior; they do not claim equivalence to the
selected Java headerless raw-byte reader or to LesionTOADS segmentation.

Two expected warnings were emitted by the synthetic Dura test because both
deliberately supplied inputs use LPS orientation and are reoriented to RAS.
No actual subject imaging data was used.

Subsequent report-only edits do not modify implementation or tests.
