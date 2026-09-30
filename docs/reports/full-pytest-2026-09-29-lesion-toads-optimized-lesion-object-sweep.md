# Full pytest evidence: P39zze OPTIMIZED_LESION OBJECT sweep

- Date: 2026-09-29
- Implementation commit tested: `312add6acbbefd7e570fa1a18c702432269bb71d`
- Command: `.venv/bin/pytest -q`
- Result: **528 passed**, 2 expected RAS orientation warnings, in 78.00 seconds.
- Focused optimized and standard OBJECT+LESION sweep suite: 8 passed.

Two committed selected-Java synthetic fixtures test `OPTIMIZED_LESION` in the
middle and last ordinary-class positions. They compare every checked interior
class and lesion float32 membership bit, including the source smoothing,
shape-gated conversion, shared normalization, and unchanged boundaries. A
separate paired Python control confirms exact alias identity with nonzero
atlas priors and distance modulation. The Java fixtures cover the
smoothing/normalization slice; prior and distance operators have separate
selected-Java references. This is not whole-engine equivalence.
