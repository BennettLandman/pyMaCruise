# Full-suite evidence: TOADS LPS source-axis compatibility

The complete repository test suite passed on implementation commit
`0981b8d153b700615450ffa64844b2404d584137`, which adds opt-in compatibility
with axial LPS source voxel order while preserving RAS image outputs.

Command:

```bash
PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result: **1,573 passed, 4 skipped, 0 failed**, with 18 expected orientation
warnings, in **708.85 seconds** (11 minutes 48 seconds). Peak sampled RSS was
**14.857 GB** under the 16 GB watchdog; the run exited normally and was not
terminated for memory.

This suite covers the Python regression suite on that commit. The source-axis
mode's case 100249 Java comparison is documented separately in the [natural
parity report](human-toads-100249-lps-source-mode-natural-java-parity-2026-10-14.md).
Neither result establishes downstream human ribbon acceptance or full
historical MultiAtlas/SPECTRE equivalence.
