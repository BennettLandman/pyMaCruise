# Full pytest evidence — 2026-09-29 LesionTOADS field data energies

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest -q` from the repository root |
| Result | **431 passed in 54.22 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `c43ba3aaf7f38582dbd087bdb971f0d91bef2ebf` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39v matches nine selected-method Java formula cases for the field-corrected
`OBJECT`, `LESION`, and `BLACKHOLE` data-term branches. The cases compare the
class, lesion, and black-hole energy accumulators for three types over three
multichannel synthetic image/field vectors with nonuniform channel weights.
The field is multiplied by each modality image before centroid residuals,
preserving the source's float32 expression order. The focused P39v tests
report **2 passed**.

This validates field application and raw data energies only. Special
Sulcal-CSF/MPRAGE and FLAIR PowerTable handling, smoothing, atlas priors,
inversion, and membership normalization remain outside this harness. See
[`tests/data/synthetic/lesion_toads_field_data_terms/README.md`](../../tests/data/synthetic/lesion_toads_field_data_terms/README.md).
