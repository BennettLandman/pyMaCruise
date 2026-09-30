# Full pytest evidence: P18bl TOADS boundary-distance update

- Date: 2026-09-29
- Implementation commit tested: `5341266f146ef305c2412f31d846b74dcff01d1d`
- Command: `.venv/bin/python -m pytest -q`
- Result: **560 passed**, 2 expected RAS orientation warnings, in 77.53 seconds.
- Focused Java comparison and validation checks: 6 passed.

Selected Java and Python outputs match every signed-byte field value for
fresh and prepopulated synthetic states at limits 0, 1, and 3. The comparison
documents the selected source's non-propagating distance behavior; it does
not establish that downstream distance-weighted membership is active or that
the full TOADS algorithm is equivalent.
