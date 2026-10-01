# Full test suite: startup thinning, LesionTOADS competition, and CLI coverage

- **Command:** `pymacruise/.venv/bin/pytest -q` (run from `/Users/landmaba/dev/mipav`)
- **Implementation commit tested:** `da4524bb725c18dc1a922fb48c056f4b14daf93a`
- **Result:** 953 passed, 0 failed, 0 skipped
- **Warnings:** 11 expected RAS reorientation warnings (generic calculator, N4, Dura topology, and workflow-preflight synthetic cases)
- **Duration:** 692.89 seconds (11:32.89)
- **Date:** 2026-09-30

This run includes the selected TOADS startup relation/fielded-membership to
first-thinning composition on a folded 25³ RAS phantom, and the selected
LesionTOADS narrow-band competition pass for both LESION and BLACKHOLE
`lesionClass` settings. The competition reference uses the selected Java
18/6 topology LUT and compares the complete insertion/pop/stale/topology/
requeue trace, final state, class-volume bits, and unchanged input state.

The suite also covers repeated standard-mode LesionTOADS continuation and
convergence, initial alignment bounds, the config-driven workflow preflight,
and standalone topology-correction and surface-guided refinement CLIs. The
refinement CLI retains `source_faithful=true` by default and tests corrected
modal-label behavior only when explicitly selected. P35 recorded zero skips,
and its optional ANTsPy numerical contract ran in this environment.

All test fixtures are committed. The fixed atlas assets are bundled algorithm
resources rather than synthetic phantoms or subject test cases; the inventory
records that distinction. This suite validates the listed components and
synthetic compositions only; it does not establish complete TOADS,
LesionTOADS, MaCRUISE, or MaCRUISEplus workflow equivalence. See the relevant
plan rows for each scope boundary.
