# P18cd full-suite evidence: default five-pass Dura wrapper

- **Implementation commit:** `1488bd4573162575ac2c20eb1b7032045190dca7`
- **Command:** `.venv/bin/pytest -q`
- **Result:** 585 passed, 2 expected RAS reorientation warnings, in 78.52 seconds.
- **Focused check:** `.venv/bin/pytest -q tests/test_toads_dura.py -k 'composes_selected_stages or source_default_five'` — 2 passed, 24 deselected, in 0.49 seconds.
- **Synthetic case:** the committed 9³ Dura topology fixture supplies selected-Java-reference WM and GM topology outputs. The test omits `max_iterations`, exercises the wrapper's default of five outer TGDM passes, and compares all level-set and final masked products to a separately assembled composition on that grid. Outputs retain RAS geometry.
- **Limit:** the complete Dura wrapper is not compared to one Java executable or harness. The topology and TGDM stages have separate selected-Java comparisons; this evidence validates their five-pass wrapper orchestration and final product assembly.

Fixture provenance is recorded in [`tests/data/synthetic/dura_topology_two_pass/README.md`](../../tests/data/synthetic/dura_topology_two_pass/README.md). The earlier one-pass wrapper composition report is [P18cb](full-pytest-2026-09-29-dura-removal-composition.md).
