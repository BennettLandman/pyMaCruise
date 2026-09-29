# Full pytest evidence — 2026-09-29 topology comparison matrix

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **394 passed in 62.95 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `5bb92b79c3a80e6fd57bebfc2b6c10505e99db20` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

The run includes five committed selected-Java topology fixtures and the first
Dura WM topology adapter. The Java and Python image results match exactly for
all 3,807 voxels. Per-case event traces match all 4,648 queue-pop coordinates
and float32 priorities, and agree on all 1,511 critical-candidate decisions.
See [the comparison report](topology-correction-java-comparison-2026-09-29.md)
for phantom definitions and harness provenance.

The run also includes the P39e LesionTOADS shape-prior `.nii.gz` loader.
P39e validates descriptor order, values, aliases, RAS normalization, and
invalid path/dimension handling; it does not include registration or
transformed-prior sampling.

Two warnings are expected: a Dura test deliberately supplies LPS-oriented
volume and paint inputs and checks that they are warned about and reoriented
to RAS. No clinical or user imaging data was used. The result validates this
synthetic test suite, not a complete MaCRUISE workflow or clinical equivalence.
