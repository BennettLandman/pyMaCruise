# P18cb full-suite evidence: composed RAS Dura removal stages

- Implementation commit tested: `dabacf43cb49cfe80847a86efbd5b4a44e72aaa5`
- Command: `.venv/bin/pytest -q`
- Result: **583 passed, 2 warnings in 86.46s**
- Warnings: expected input-orientation warnings from synthetic Dura topology
  tests; the LPS test volume and paint mask are reoriented to RAS.
- Focused result: Dura module — **25 passed**.

`run_dura_removal` composes the selected first WM topology call, float32
WM-plus-raw-GM union and second topology call, initial signed-distance fields,
one configurable TGDM outer pass, and final threshold/output masking. The two
intermediate topology volumes match the selected Java harness; the assembled
products match the selected synthetic stage sequence and remain on the T1 RAS
grid. This validates wrapper composition on a synthetic 9³ case. Full default
five-pass and folded-case wrapper acceptance, plus complete TOADS, remain
open.
