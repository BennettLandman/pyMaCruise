# Full test suite: continuation and source-evidence batch

## Run record

- **Implementation commit:** `c3ad19daefe7cd0987b1158fa48b43e514a405c2`
- **Command:** `.venv/bin/pytest -q`
- **Result:** **922 passed, 4 warnings, 520.59 seconds (8:40)**
- **Failures/skips:** none reported

The four warnings are expected RAS conversion notices: two are emitted when
the calculator test constructs its expected output from LAS synthetic inputs,
and two are emitted when the Dura topology test normalizes LPS topology and
paint-mask inputs to RAS.

This complete run includes the new selected-source continuation slices,
synthetic Java references, and provenance assertions. It does not establish
complete TOADS/LesionTOADS workflow equivalence, real-data readiness, or P37
acceptance.

## Focused evidence in the full run

- P18cyk: the four-class generic continuation/thinning chain and selected
  11-class next-thinning fixture both compare source-ordered Java references.
- P39zzzzzzal: refreshed centroids, class/auxiliary standard deviations, and
  channel weights compare bit-for-bit with a Java harness on the exact 17³,
  13-class post-thinning state.
- P39zzzzzc, P39zzzzy, P39zzzzzd, P39zzzzze, P39zzzzzaa, P39zzzzzae,
  P39zzzzzaf, P39zzzzzzag, and P39zzzzzzah: source/harness/input hashes,
  wrapper order, source traversal, signed-byte behavior, and full synthetic
  output/reference assertions are exercised by their row-specific tests.

