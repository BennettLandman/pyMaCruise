# Full pytest evidence — 2026-09-29 LesionTOADS registration pyramid

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **410 passed in 54.14 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `a194c16c3fe97aaacdb21797d8d85f0c79d19912` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39o composes selected source block-mean downsampling with fixed-scale
LesionTOADS registration optimization. On an odd 33³ two-class synthetic
case, the Java harness and Python match coarse dimensions, scale order,
transform handoff, every line-search trial, stage transform/energy/lambda,
and final transform for levels 1 through 4. The selected method omits scale 1;
levels 1 and 4 produce no stage and scales 8→4→2 respectively. Some large
scale-8 trial transforms yield a Java NaN energy at zero prior norm; JSON
`null` in the fixture represents that NaN, and both paths reject those trials
before later accepting a valid step. Provenance is in
`tests/data/synthetic/lesion_registration_pyramid/README.md`.

This is selected-formula pyramid registration evidence only, not a full
`LesionToadAtlas` run, lesion segmentation, or plus-workflow equivalence. The
full suite used synthetic cases; the two RAS warnings are expected orientation
checks from Dura topology tests.
