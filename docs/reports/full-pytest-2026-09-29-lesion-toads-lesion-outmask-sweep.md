# Full pytest evidence: P39zy LESION+OUTMASK volume composition

- Date: 2026-09-29
- Implementation commit tested: `c712c2ab983aa406d1983619b7c08ea4e446ba19`
- Command: `.venv/bin/pytest -q`
- Result: **509 passed**, 2 expected RAS orientation warnings, in 76.82 seconds.
- Focused combined/adjacent LesionTOADS volume suite: 11 passed.

The committed synthetic `5×4×6` selected-Java formula fixture composes the
ordinary class, LESION, and OUTMASK paths in one source-ordered membership
pass. It includes class and lesion atlas similarities, six-face smoothing,
mixed direct/prior relations, a zero-shape lesion voxel, zero/unit auxiliary
weights, and shared denominator normalization. Every interior class, lesion,
and outlier float32 membership bit matches the Java harness; outer planes are
unchanged, inputs are not mutated, and zero atlas-prior coefficient matches
the combined no-prior control. This is a selected volume formula validation,
not a full LesionTOADS-equivalence claim.
