# Full synthetic pytest: TOADS atlas through post-refresh registration

**Tested implementation commit:** `ee350877e61f97adb050d185bfb34f7bfd067ade`  
**Date:** 2026-10-02  
**Command:** `./.venv/bin/pytest -q`

Result: **1,049 passed, 0 failed, 0 skipped, 11 expected RAS orientation
warnings in 1,682.82 seconds**.

This major synthetic integration checkpoint includes atlas descriptor and RAS
asset loading, initial rigid startup, post-registration native-grid recropping
and topology refresh, and exact-membership initialization followed by the
selected single-scale rigid registration update. The checkpoint does not
establish complete TOADS or MaCRUISE equivalence. It does not use the supplied
human MRI run.

The selected TOADS source-line representation estimate was re-evaluated with
this checkpoint and remains **3,400–4,800 of 13,651 lines (about 25–35%)**.
The interval stays unchanged because the latest adapter composes operation
families already represented by existing validated source slices; this is a
coarse method-level estimate, not runtime line coverage. See the [source-line
accounting note](toads-source-line-accounting-2026-10-01.md) and the current
estimate at the top of [plan.md](../../plan.md).
