# Full pytest evidence — 2026-09-29 scaled LesionTOADS optimizer loop

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **409 passed in 53.07 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `76e6ec46a3940a3e5a584e725a9c751962e9e32e` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39n implements the selected fixed-scale outer loop using P39m steps. On the
committed deterministic 8³ synthetic case, Java and Python match every
iteration's trial transforms and energies, accepted transform, energy,
lambda, relative change, and rejection count. The three harness controls stop
by source lambda floor after six iterations, by a three-iteration cap, and by
a minimum relative-change threshold after one iteration. Fixture provenance
and regeneration instructions are in
`tests/data/synthetic/lesion_scaled_registration_optimizer/README.md`.

This validates one scale's iteration/stopping loop only. It does not construct
or transition the pyramid, optimize transformed bounds, perform deformable
registration, or segment lesions. The full suite used synthetic cases; the
two RAS warnings are expected orientation checks from Dura topology tests.
