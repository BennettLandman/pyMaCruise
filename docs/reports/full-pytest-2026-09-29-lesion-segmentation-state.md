# Full pytest evidence — 2026-09-29 LesionTOADS segmentation state

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **416 passed in 54.01 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `c2b14be265412e7d011c02402a4f0f3bf4a78131` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39p adds constructor-time state initialization for LesionTOADS. The compact
synthetic RAS fixture and selected-constructor Java formula harness check
class-index segmentation, skeleton and relation maps; class volumes and
source-order float32 image centroids; zero membership/availability; empty
class NaN centroids; lesion/no-lesion class similarity defaults; configured
distance sentinels; and other initial model arrays. The focused P39p tests
report **6 passed**.

This is constructor-state formula evidence only. The harness does not
instantiate the entire VandyCRUISE class; no lesion membership evolution,
topology propagation, segmentation, or plus-workflow equivalence is claimed.
Fixture provenance and image checksums are in
`tests/data/synthetic/lesion_segmentation_state/README.md`.
