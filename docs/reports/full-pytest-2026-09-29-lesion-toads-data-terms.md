# Full pytest evidence — 2026-09-29 LesionTOADS unfielded data energies

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest -q` from the repository root |
| Result | **429 passed in 54.25 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `c69df55991a468279f1e31b213cb64783f910f24` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39u matches nine selected-method Java formula cases for the unfielded
`OBJECT`, `LESION`, and `BLACKHOLE` data-term branches. The cases compare the
class, lesion, and black-hole energy accumulators for three types over three
multichannel synthetic samples with nonuniform channel weights. The Python
operation order preserves Java float32 accumulation. The focused P39u tests
report **2 passed**.

This validates raw intensity energies only. Bias-field multiplication,
Sulcal-CSF/MPRAGE or FLAIR PowerTable handling, smoothing, atlas priors,
inversion, and membership normalization are outside the harness. See
[`tests/data/synthetic/lesion_toads_data_terms/README.md`](../../tests/data/synthetic/lesion_toads_data_terms/README.md).
