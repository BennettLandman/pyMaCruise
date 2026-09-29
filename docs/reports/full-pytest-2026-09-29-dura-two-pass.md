# Full pytest evidence — 2026-09-29 Dura two-pass topology

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **395 passed in 64.08 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `b2ac696ee453880eabcc1d85dfcf22d6c6e4d9c0` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

The run includes the P18bf synthetic comparison for the first Dura WM
topology output, the source-order float32 union with raw GM, and the second GM
topology output. All three match the selected Java algorithm-pair harness
voxelwise exactly on the 9×9×9 fixture, with float32 dtype and identical RAS
geometry. The harness executes the two selected topology algorithms and union
as a test oracle following `MedicAlgorithmRemoveDura09`; it does not instantiate
the full JIST wrapper. Grid mismatch rejection is covered separately.

The run also includes five exact P25 topology phantoms and all previous tests.
Two expected warnings come from a synthetic LPS Dura test and confirm that its
volume and paint mask are reoriented to RAS. No clinical or user imaging data
was used. This does not establish complete Dura workflow or TOADS equivalence.
