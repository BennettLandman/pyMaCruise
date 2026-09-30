# Full pytest evidence — atlas and LesionTOADS branch follow-up

**Run date:** 2026-09-30  
**Command:** `.venv/bin/pytest -q`  
**Implementation commit:** `3accee53bb1ab1c8787e5af8791348cd77a373dc`  
**Result:** **796 passed, 2 warnings, 154.78 seconds**

The two warnings are expected RAS-normalization warnings from
`test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`:
the synthetic topology volume and paint mask are LPS-oriented and are
explicitly reoriented to RAS. No test failed.

This full synthetic suite validates the committed implementation slices:

- **P18czj:** continuing default TOADS state into the next thinning and
  fixed-statistics membership pass. Its 25³ folded phantom produced 8,962
  queue pops and 978 source-matched segmentation updates; next-pass energies,
  memberships, denominator, and objective matched the selected Java reference.
- **P18czk:** the bundled 13-structure TOADS2012 exact-membership dispatch,
  including all 16 field/PSF/smoothing/prior reference combinations.
- **P39zzzzu:** the LesionTOADS post-registration stop/continue branch using
  P39zzzzs/P39zzzzt state and the selected Java wrapper control/order reference.
- The earlier completed slice rows also remain covered by this repository-wide
  suite; their first full-run evidence is cited individually in `plan.md`.

The run is on committed implementation snapshot `3accee5`. It did not include
the in-progress, untracked P39zzzzv relation-aware next-iteration draft; that row
remains Started. All images exercised here are synthetic. This run does not
claim complete TOADS or LesionTOADS workflows or real-data acceptance.
