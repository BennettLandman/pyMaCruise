# Full pytest evidence: P39zzd OPTIMIZED_LESION+OUTMASK volume composition

- Date: 2026-09-29
- Implementation commit tested: `6fecb9cdf48573d6b61598c009c325a38249a699`
- Command: `.venv/bin/pytest -q`
- Result: **524 passed**, 2 expected RAS orientation warnings, in 84.97 seconds.
- Focused LesionTOADS auxiliary sweep suite: 15 passed.

Two committed selected-Java synthetic references replace LESION with
OPTIMIZED_LESION and exercise both positions relative to OUTMASK. They compare
every interior class, lesion, and outlier float32 membership bit, including
smoothing and atlas priors, shared normalization, zero shape support, and
boundary preservation. The OUTMASK prior fixture preserves the Java expression
order `shapeFactor(shape[k]) * outlierMembership * numO`. This validates the
selected combined formula path, not the complete iterative LesionTOADS
workflow.
