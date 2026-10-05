# Full-suite evidence: restartable plus LesionTOADS loop workflow

- **Tested implementation commit:** `c74a739107fe814d9fdec1f632b3fe1af5d7fbb2`
- **Full command:** `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- **Result:** 1,188 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 982.36 seconds.
- **Memory:** Peak sampled RSS 15.112 GB; the 24 GB watchdog did not terminate the run.
- **Focused P41o chain:** `PYTHONPATH=src ./.venv/bin/pytest -q tests/test_plus_lesion_registered_branch.py` — 8 passed in 72.12 seconds.
- **CLI help:** `PYTHONPATH=src ./.venv/bin/python -m pymacruise.cli plus-lesion-loop-run --help` lists the recrop/atlas inputs, numeric settings, bounded pass, resume, source-faithful default, and explicit finite-mode switch.

P41o adds one command from a completed P41h registered recrop through P41i initial membership/thinning, P41k registered control using P41i's persisted state, and a bounded P41m loop checkpoint. A top-level JSON manifest records running/complete stage states and hashes for all child manifests; each child's own checkpoint code verifies output hashes on resume. The synthetic integration test checks RAS child images, iteration-2 persistence, valid resume without recomputation, and repair of a corrupt loop output while reusing verified predecessors.

This is synthetic selected LesionTOADS continuation through one bounded loop segment. It does not produce terminal LesionTOADS outputs, complete MaCRUISEplus, establish composite Java equivalence, or use human imaging data. No new selected Java methods were implemented, so the coarse TOADS source-operation estimate remains 3,464–4,864 of 13,651 lines (25–36%), with a zero line delta for this integration milestone.
