# Full synthetic pytest: TOADS initial membership through first thinning

**Tested implementation commit:** `b6b231d8f3eeed914187c577b6423b6e1fb1681a`  
**Date:** 2026-10-02  
**Command:** `./.venv/bin/pytest -q`

Result: **1,054 passed, 0 failed, 0 skipped, 11 expected RAS orientation
warnings in 1,826.52 seconds (30:26)**.

This checkpoint includes the committed startup composition through initial
class deviations, image-wide variance weights, centroid similarities, one
exact membership/energy pass, and first thinning transition, alongside the
previous synthetic repository coverage. The composition still supplies
transformed atlas shape samples and other upstream values explicitly; this
full test run is not evidence of complete TOADS or MaCRUISE equivalence and
does not use the supplied human MRI run.

The TOADS source-line representation estimate was re-reviewed with this suite
and remains **3,400–4,800 of 13,651 selected source lines (about 25–35%)**. The
new adapters connect already represented relation, deviation, image-variance,
similarity, membership, and thinning operations. They do not add a distinct
source method family to the coarse method-level numerator. This is not runtime
Java line coverage; see the [source-line accounting note](toads-source-line-accounting-2026-10-01.md)
and [current estimate](../../plan.md).
