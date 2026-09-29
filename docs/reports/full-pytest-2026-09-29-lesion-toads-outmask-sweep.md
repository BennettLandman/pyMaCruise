# Full pytest evidence: P39zw OUTMASK volume sweep

- Date: 2026-09-29
- Implementation commit tested: `6895cf32d9af1ba989af2bd069997ca4bc62092f`
- Command: `.venv/bin/pytest -q`
- Result: **504 passed**, 2 expected RAS orientation warnings, in 78.06 seconds.
- Focused OUTMASK/normalization/smoothing suite: 5 passed.

The committed synthetic `5×4×6` Java formula fixture exercises source-ordered
in-place smoothing, one OUTMASK class, mixed direct/prior relations, zero,
interior, and unit outlier weights, shape gating, outlier-map updates, and
shared denominator normalization. Every interior class-membership and
outlier-membership float32 bit matches the selected Java harness; boundaries
remain unchanged. This validates the selected OBJECT+OUTMASK branch only; it
does not claim whole LesionTOADS equivalence or cover atlas similarity priors.
