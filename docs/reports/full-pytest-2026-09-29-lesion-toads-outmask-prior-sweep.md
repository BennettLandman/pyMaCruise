# Full pytest evidence: P39zx OUTMASK atlas-prior volume sweep

- Date: 2026-09-29
- Implementation commit tested: `ab80898125b2640f6f4bb118e4a6981448285ea5`
- Command: `.venv/bin/pytest -q`
- Result: **506 passed**, 2 expected RAS orientation warnings, in 77.67 seconds.
- Focused OUTMASK/normalization/smoothing suite: 7 passed.

The committed synthetic `5×4×6` Java formula fixture adds asymmetric atlas
shape priors and off-diagonal class similarities to the OBJECT+OUTMASK volume
sweep. Every interior class and outlier-membership float32 bit matches the
selected Java calculation, including the shared class-prior contribution to
`numO`. A zero atlas-prior coefficient exactly reproduces P39zw, and boundary
voxels remain unchanged. This validates the selected formula branch only; it
does not claim whole LesionTOADS equivalence.
