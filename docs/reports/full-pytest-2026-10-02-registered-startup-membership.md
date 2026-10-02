# Full synthetic test suite — registered TOADS startup membership

- **Tested implementation commit:** `dcb707fec0425c3e28fe7b20bf2b2c5822a061be`
- **Command:** `./.venv/bin/pytest -q`
- **Result:** **1,060 passed, 0 failed, 0 skipped; 11 expected RAS-orientation warnings**
- **Elapsed time:** **1,862.37 seconds (31:02.37)**
- **Date:** 2026-10-02

This suite includes the P18cyat adapter and the new synthetic integration
which builds a selected 11-class atlas from committed-in-test `.nii.gz`
assets, runs `register_toads_from_atlas_assets`, and carries that actual
`ToadsAtlasAlignmentRun` into refreshed statistics and an exact membership
pass. The suite passes without MATLAB, MIPAV, JIST, FreeSurfer, or real subject
data.

The integration compares the composed membership array, denominator, and
objective against a direct call to the existing exact-membership kernel using
the run's registration-derived arrays. The component kernels have separate
selected-Java reference coverage. This remains a synthetic adapter comparison,
not a single combined Java wrapper oracle or a complete TOADS iteration.

This checkpoint added a tested startup-to-membership integration adapter but
no new standalone Java algorithm method family to the source-line accounting.
The source-line estimate was therefore re-reviewed and remains
**3,400–4,800 of 13,651 selected Java lines (about 25–35%)**. It remains a
coarse represented-operation estimate, not runtime line coverage.
