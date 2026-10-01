# Full test suite: TOADS and LesionTOADS continuations, N4 CLI

- **Command:** `.venv/bin/pytest -q`
- **Implementation commit:** `07d4098db791d03db5b0fd90fd1d6794123897bd`
- **Result:** 933 passed, 0 failed
- **Warnings:** 5 expected RAS reorientation warnings (calculator inputs, N4 input, and Dura topology image/mask inputs)
- **Duration:** 578.92 seconds (9:38.92)
- **Date:** 2026-09-30

This run covers the source-matched 11-class TOADS continuation through the
next-loop thinning, centroid/variance refresh, and two-channel exact-membership
update; the no-atlas and atlas-aware TOADS approximate-membership references and
initial registration handoff; the LesionTOADS post-thinning `GMflag=true`
membership pass with Java-compatible zero-variance NaNs and the following
float32 energy-ratio calculation; and the standalone N4 CLI with its optional
ANTsPy dependency path. The synthetic fixtures are under
`tests/data/synthetic/` and contain no subject images.

The ratio harness records `Numerics.abs(float)` behavior directly, including
the selected Java NaN result `0xffc00000` and negative-zero ratio `0x80000000`.
The initial TOADS registration comparison remains bounded to the documented
states and measured float differences in P18czq; this suite does not establish
full TOADS, LesionTOADS, MaCRUISE, or MaCRUISEplus workflow equivalence. See
the corresponding plan rows for those scope boundaries.
