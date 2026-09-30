# Full pytest evidence: P39zzr LesionTOADS channel weights

- Run: `.venv/bin/pytest -q`
- Result: **614 passed, 2 warnings, in 90.66s**.
- Warnings: the two expected RAS-normalization warnings from tests that deliberately supply LPS-oriented NIfTI images.
- Tested implementation commit: `79ad6140610123725a1d758e4602aac8592a9b19`.
- Focused selected-Java channel-weight comparisons: `tests/test_lesion_channel_weights.py` (**2 passed**).
- Synthetic input and reference details: [channel-weight phantom](../../tests/data/synthetic/lesion_channel_weights/README.md).

The Python implementation and selected-method Java excerpt match raw float32
channel-weight bits for the diagonal-covariance mode with and without field
correction. The synthetic volume exercises OBJECT, OUTLIER, LESION, and
BLACKHOLE contributions, while MASK/OUTMASK are present as source-ignored
types. It also covers the source's linear lesion-related residual terms, the
uniform fallback for other modes, and the nonpositive-energy path where each
channel uses `INF=1e10` before normalization to equal weights.

This validates channel weights for fixed inputs only; it does not validate
the full class centroid and membership-estimation pipeline.
