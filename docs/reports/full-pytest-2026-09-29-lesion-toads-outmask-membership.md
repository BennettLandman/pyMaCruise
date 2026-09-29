# Full test suite: LesionTOADS OUTMASK membership

- Date: 2026-09-29
- Tested implementation commit: `94a919564650fddd882e3cb30ba3dc6763e1a690`
- Command: `.venv/bin/pytest -q`
- Result: **453 passed in 57.68 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes twelve selected-Java scalar formula cases for the
`OUTMASK` mixture. Direct/prior class output and updated outlier membership
bits match across shape factors below, at, and above saturation, outlier
weights from zero to one, and float32-sensitive values. The fixture starts
from already inverted class/outmask terms; denominator normalization is not
included in this test slice.
