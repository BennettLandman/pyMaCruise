# Full pytest evidence — 2026-09-29 LesionTOADS membership conversion

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest -q` from the repository root |
| Result | **438 passed in 54.59 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `6c92a3db7da65c96b2f01c0f06e791cd479b9813` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39y matches eleven selected-Java formula cases for class energy inversion
and direct/prior membership conversion. Cases test below, exactly at, and
just above the strict `ZERO` cutoff; zero/negative energies mapping to the
finite `INF` sentinel; and `shapeFactor = Numerics.min(1, 2*shape)` below,
at, and above saturation. Three focused tests pass.

This validates one class conversion only. It does not calculate lesion
auxiliary membership, distance/relation modulation, denominator accumulation,
or voxelwise normalization. Fixture details are in
[`tests/data/synthetic/lesion_toads_membership_conversion/README.md`](../../tests/data/synthetic/lesion_toads_membership_conversion/README.md).
