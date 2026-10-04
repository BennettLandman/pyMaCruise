# Full pytest after active LesionTOADS 13-class raw energy dispatch

Command:

```text
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q
```

Tested commit: `d60c99777b2a776b95a7898b4120918559c935b0`.

Result: **1,146 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
804.55 seconds (13:24)**. The RSS watchdog exited 0, sampled peak RSS of
**15.764 GB**, and did not terminate the child process (24 GB cap).

The suite used synthetic data and did not execute the supplied human MRI.
P39zzzzzzbd adds source-bit-matched raw data terms for the active classes of
the bundled 13-class LesionTOADS descriptor. It does not complete membership
updates or the LesionTOADS workflow. No new source lines were counted toward
the selected base-TOADS denominator (3,464–4,864 of 13,651); P39 is tracked
separately. This is test-regression plus bounded source-formula evidence, not
end-to-end workflow equivalence or real-data validation.
