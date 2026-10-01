# Full synthetic pytest run — 2026-10-01

Implementation commit under test: `6697e0c8e12d28b89af43680665a385cee624133`.

Command, run from the `pymacruise` repository root:

```text
./.venv/bin/pytest -q
```

Result: **976 passed, 0 failed, 0 skipped, 11 expected RAS-reorientation warnings** in **699.12 seconds (11:39)**.

The 11 warnings are expected checks that non-RAS inputs are reoriented: two image-calculator inputs, one N4 input, two Dura topology inputs, and six workflow-preflight inputs. No test warning represents a numerical mismatch or failed test.

This suite includes the selected P18 registration-through-outer-loop composition, P39's bounded 13-class source-generated-NaN thinning pass, P36a's reference-corpus index, and the isolated-wheel P42a2 preflight. P39's NaN mode is intentionally limited to source events skipped by Java's volume gate; it does not exercise a NaN-dependent topology edit. A post-run code review found that an unsupported active-NaN event pops the caller-owned queue before the guard raises; P39as remains Started until that error-path behavior is fixed and revalidated.

No actual subject imaging data was used. Bundled atlas/template images remain committed fixed algorithm resources, not synthetic test cases. The suite does not establish complete selected TOADS, LesionTOADS, base MaCRUISE, or plus workflow equivalence; those remain open in [plan.md](../../plan.md).
