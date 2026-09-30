# P18bz full-suite evidence: TOADS OBJECT/OUTLIER atlas-prior composition

- Implementation commit tested: `34137d63ee76781cbcc62fa2287e500fff768d8c`
- Command: `.venv/bin/pytest -q`
- Result: **581 passed, 2 warnings in 77.14s**
- Warnings: expected input-orientation warnings from synthetic Dura topology
  tests; the LPS test volume and paint mask are reoriented to RAS.
- Focused result: selected Java OBJECT/OUTLIER sweep comparison with and
  without atlas inverse-prior energy — **1 passed**.
- Java reference generator checks the selected
  `MultipleToadSegmentation.java` SHA-256 before compiling the harness.

The fixture matches every inverse-prior and half-weighted objective-energy
bit, then every final class membership, scalar outlier-map, and total-objective
bit. The compared cases are synthetic and cover the selected membership pass;
they do not establish topology admissibility or whole-TOADS equivalence.
