# Full pytest evidence — 2026-09-29 scaled LesionTOADS registration step

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **408 passed in 52.49 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `f2216f08ee75a861eba86e6c4e7d5f6103a64c95` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39m composes fixed-transform coefficient calculation with the existing
source-checked diagonal gradient update and strict-improvement line search.
Each candidate recomputes the scaled LesionTOADS volume energy from the same
synthetic atlas maps and memberships. At scales 1 and 2, the selected-Java
formula harness and Python implementation match every trial transform and
energy, plus the final transform, energy, lambda, relative change, accepted
flag, and rejection count. The scale-1 case accepts immediately; the scale-2
case backtracks eight times before acceptance. Fixture provenance and
regeneration steps are in
`tests/data/synthetic/lesion_scaled_registration_step/README.md`.

This is one fixed-scale registration update, not an outer iteration loop,
convergence check across iterations, deformable registration, or lesion
segmentation. The full suite used synthetic cases; the two RAS warnings are
expected orientation checks from Dura topology tests.
