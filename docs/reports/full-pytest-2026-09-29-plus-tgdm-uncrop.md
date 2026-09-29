# Full test suite: MaCRUISEplus TGDM output uncropping

- Date: 2026-09-29
- Tested implementation commit: `699c25aa550c6cd9315bfc44d392a345aeaf0c7e`
- Command: `.venv/bin/pytest -q`
- Result: **482 passed in 78.15 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

`uncrop_plus_tgdm_outputs` restores the six scalar TGDM products through the
existing JIST-compatible crop helper using the GM-derived `CropParameters`.
On the synthetic 9³ reference with known offset and anisotropic RAS affine, all
six arrays match direct `uncrop_image(..., legacy_exclusive_max=True)` and
return to the full reference shape and affine. The test checks the selected
source behavior that omits the maximum crop plane and rejects an image on a
different cropped grid.

Focused tests: 6 passed in 0.53 seconds, including existing crop semantics.
This restores TGDM volume outputs; mesh vertex offset conversion and complete
plus workflow orchestration remain separate work.
