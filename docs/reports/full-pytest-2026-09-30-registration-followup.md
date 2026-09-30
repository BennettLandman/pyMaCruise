# Full pytest evidence — registration handoff follow-up

**Run date:** 2026-09-30  
**Command:** `.venv/bin/pytest -q`  
**Implementation commit:** `3e1d8f1d3d16f2851c0b0297654c7d154467987d`  
**Result:** **785 passed, 2 warnings, 131.02 seconds**

The two warnings are expected RAS-normalization warnings from
`test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`:
the synthetic topology volume and paint mask are stored in LPS orientation and
are explicitly reoriented to RAS. No test failed.

This full synthetic suite includes the source-pinned method and workflow-slice
tests for the following roadmap rows:

- **P18cx:** corrected OUTLIER auxiliary smoothing against the selected Java
  dangling-`else` behavior.
- **P18czd–P18czg:** LESION, OUTMASK, default 11-class dispatch, and OUTLIER
  membership compositions.
- **P18czh–P18czi:** default membership-to-registration handoff and immediate
  post-registration stop/grow/classification-or-refresh branch.
- **P39zzzzq–P39zzzzt:** dormant covariance statistics, active statistics
  refresh, statistics-to-membership composition, and membership-to-registration
  handoff.

The source-slice reports and committed Java harness references define each
row's numerical scope. This run establishes that the full repository test suite
passes at the implementation snapshot above; it does not claim that the full
TOADS or LesionTOADS workflows, MaCRUISEplus, or real-data acceptance are
complete. All imaging inputs exercised by this run are synthetic.
