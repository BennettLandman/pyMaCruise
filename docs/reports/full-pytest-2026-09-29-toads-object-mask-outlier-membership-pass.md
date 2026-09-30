# P18ca full-suite evidence: TOADS OBJECT/MASK/OUTLIER membership pass

- Implementation commit tested: `22f1bc4d8267af64f0cd204dc5bf6a1d0c8c28f7`
- Command: `.venv/bin/pytest -q`
- Result: **582 passed, 2 warnings in 77.11s**
- Warnings: expected input-orientation warnings from synthetic Dura topology
  tests; the LPS test volume and paint mask are reoriented to RAS.
- Focused result: selected Java OBJECT/MASK/OUTLIER composed-pass comparison —
  **1 passed**.
- Java reference generator checks the selected
  `MultipleToadSegmentation.java` SHA-256 before compiling the harness.

On the committed 4 x 4 x 4 synthetic case, Python and the source-formula
harness match all class-membership, scalar outlier-map, and objective
float32 bits after dispatching OBJECT residual, MASK mismatch, and paired
OUTLIER energies through the smoothing, atlas-prior, and relation-weighted
update. This is selected membership-pass evidence; it does not validate
topology gates or complete TOADS evolution.
