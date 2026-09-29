# Full pytest evidence: P39zzc OPTIMIZED_LESION volume composition

- Date: 2026-09-29
- Implementation commit tested: `6c0835ede94ec411b0d22a82cbcc2c0045d6754d`
- Command: `.venv/bin/pytest -q`
- Result: **521 passed**, 2 expected RAS orientation warnings, in 84.10 seconds.
- Focused three-auxiliary and related sweep suite: 16 passed.

Two committed selected-Java synthetic references replace LESION with
OPTIMIZED_LESION and test both positions relative to BLACKHOLE. They compare
every interior class, lesion, black-hole, and outlier float32 membership bit,
including smoothing and prior terms, distance modulation, OUTMASK blending,
shared normalization, and boundary preservation. The harness follows the
selected Java class branches; this validates the combined formula path, not
the complete iterative LesionTOADS workflow.
