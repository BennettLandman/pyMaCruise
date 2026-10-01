# Full synthetic pytest run — 2026-10-01

Implementation commit under test: `163a5e74e690c629f62ccb4e55fea6c32bff72d9`.

Command, run from the `pymacruise` repository root:

```text
./.venv/bin/pytest -q
```

Result: **973 passed, 0 failed, 0 skipped, 11 expected RAS-reorientation warnings** in **717.81 seconds (11:57)**.

The 11 warnings are expected checks that non-RAS inputs are reoriented: two image-calculator inputs, one N4 input, two Dura topology inputs, and six workflow-preflight inputs. No test warning represents a numerical mismatch or failed test.

This complete suite includes the selected VandyCRUISE P18cyp registration trace, the P39ar terminal/growth references, plus boundary and atlas-selection diagnostics, active Dura and TGDM CLI tests, and all previously committed synthetic contracts. The exact row scope and remaining integration limits are recorded in [plan.md](../../plan.md). The refreshed P35 synthetic-corpus/static-scan inventory for this same implementation commit is [here](p35-synthetic-readiness-inventory-2026-10-01.md).

No actual subject imaging data was used. Bundled atlas/template images remain committed fixed algorithm resources, not synthetic test cases. This component suite does not establish the full selected TOADS, LesionTOADS, base MaCRUISE, or plus workflow equivalence gates; those remain open.
