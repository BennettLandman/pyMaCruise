# Full pytest evidence — 2026-09-29 LesionTOADS Sulcal-CSF energy

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest -q` from the repository root |
| Result | **435 passed in 53.56 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `211d6286498a9065f92cb407f4e69dd621e692b7` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39x matches nine selected Java formula cases for the ordinary-object
Sulcal-CSF data-energy branch. They cover MPRAGE PowerTable and direct T1
residuals, field-corrected forms, absent-T1 generic multichannel fallback,
and exact structure-name gating. Class-energy float32 bits match exactly.
Two focused tests pass, including invalid T1 index and missing PowerTable
exponent rejection.

This validates only this selected per-class data term. Smoothing, atlas
priors, energy inversion, normalization, and the full membership iteration
remain open. Fixture details are in
[`tests/data/synthetic/lesion_toads_sulcal_csf/README.md`](../../tests/data/synthetic/lesion_toads_sulcal_csf/README.md).
