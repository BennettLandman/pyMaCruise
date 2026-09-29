# Full pytest evidence: P39zz combined LESION distance routing

- Date: 2026-09-29
- Implementation commit tested: `afed23eca5c3daec529339c7a00ec893a7dc8e21`
- Command: `.venv/bin/pytest -q`
- Result: **512 passed**, 2 expected RAS orientation warnings, in 77.36 seconds.
- Focused combined/adjacent LesionTOADS volume suite: 10 passed.

Two committed synthetic Java references cover direct and prior LESION
relations in the combined OBJECT+LESION+OUTMASK sweep. They exercise
ventricle precedence, GM fallback for the `-1` ventricle sentinel,
zero-distance `spread`, independent brainstem modulation, mixed outlier
weights, class/lesion atlas priors, and shared normalization. Every interior
class, lesion, and outlier float32 membership bit matches; borders remain
unchanged. Zero distance maxima reproduce P39zy exactly. This is selected
distance routing evidence, not full LesionTOADS equivalence.
