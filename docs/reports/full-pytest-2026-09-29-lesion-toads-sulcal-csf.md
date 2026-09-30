# Full pytest evidence — 2026-09-29 LesionTOADS Sulcal-CSF energy

> **Superseded for P39x branch correctness.** This report records the full-suite
> run on the earlier nine-case implementation and does not validate the later
> source correction. The corrected branch now has a 10-case Java harness and
> focused regression evidence (10 passed in 0.75 seconds); a full-suite rerun
> after that correction is pending. The original run details below remain as
> historical evidence for commit `211d6286498a9065f92cb407f4e69dd621e692b7`.

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

The run above predates the corrected Sulcal-CSF branch composition. The
selected unfielded Java path adds the special T1 term and then the generic
weighted all-channel term because the apparent `else` is commented out. The
fielded path uses an active `else` and therefore keeps the T1 term exclusive;
no T1 channel selects only the generic term. The original nine-case fixture
did not test the unfielded special-plus-generic behavior and is superseded for
that claim.

After the source correction, the regenerated 10-case harness matches the
corrected Python helper bit-for-bit, including a rounding-sensitive check of
special-term-first accumulation. Focused regression across Sulcal-CSF,
unfielded and fielded data terms, FLAIR weighting, and membership composition:
**10 passed in 0.75 seconds**. `git diff --check` passed. This focused result
does not replace the pending complete-suite rerun.

This validates only this selected per-class data term. Smoothing, atlas
priors, energy inversion, normalization, and the full membership iteration
remain open. Fixture details are in
[`tests/data/synthetic/lesion_toads_sulcal_csf/README.md`](../../tests/data/synthetic/lesion_toads_sulcal_csf/README.md).
