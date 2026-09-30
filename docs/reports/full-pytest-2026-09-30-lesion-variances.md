# Full pytest evidence: P39zzq LesionTOADS variances

- Run: `.venv/bin/pytest -q`
- Result: **612 passed, 2 warnings, in 78.04s**.
- Warnings: the two expected RAS-normalization warnings from tests that deliberately supply LPS-oriented NIfTI images.
- Tested implementation commit: `5f1948658030ec17407d1d53e3875d809b57cfda`.
- Focused selected-Java variance comparisons: `tests/test_lesion_variances.py` (**2 passed**).
- Synthetic input and reference details: [variance phantom](../../tests/data/synthetic/lesion_variances/README.md).

Python and the selected-method Java excerpt match raw float32 bits for
per-class, lesion, and black-hole standard deviations with and without field
correction. The fixture covers ordinary, LESION, and BLACKHOLE classes, an
empty class, lesion/black-hole winner branches, equality ties, and zero
auxiliary weights. The Java and Python references preserve the selected
source's black-hole winner quirk: the black-hole residual is weighted by
lesion membership in that branch.

This validates the selected variance update for fixed inputs. It does not
validate centroid-category assignment, channel weights, membership
estimation, or full LesionTOADS iteration.
