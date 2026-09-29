# Full pytest evidence — 2026-09-29 LesionTOADS FLAIR lesion weighting

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest -q` from the repository root |
| Result | **433 passed in 54.03 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `7ee767ca6c2c498115561ba5f57721e2ef3047e4` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39w matches six selected Java lesion-likelihood cases using the actual
source `PowerTable`: three synthetic multichannel samples, each with and
without field correction. The FLAIR-specific squared residual includes both
cached-range and out-of-range direct-power inputs. Class-energy and
float32-narrowed lesion-energy bits match exactly (maximum absolute error
0.0). Two focused tests pass, including absent-FLAIR fallback and required
power-exponent checks.

This validates only the selected `LESION`/`useLesionWeight` data-term branch.
It does not include Sulcal-CSF/MPRAGE weighting, smoothing, atlas priors,
inversion, or normalized memberships. Fixture provenance is in
[`tests/data/synthetic/lesion_toads_flair_weight/README.md`](../../tests/data/synthetic/lesion_toads_flair_weight/README.md).
