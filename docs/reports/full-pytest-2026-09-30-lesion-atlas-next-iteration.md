# Full pytest evidence — LesionTOADS atlas and continuing iteration

**Run date:** 2026-09-30  
**Command:** `.venv/bin/pytest -q`  
**Implementation commit:** `d84f2039fe45e6c90cd040a559cc09029a379b7a`  
**Result:** **810 passed, 2 warnings, 216.16 seconds**

The two warnings are expected RAS-normalization warnings from
`test_toads_dura.py::test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`: the
synthetic topology volume and paint mask are LPS-oriented and are explicitly
reoriented to RAS. No test failed.

This complete synthetic suite covers the following newly completed slices:

- **P18czl:** the selected 13-class LesionTOADS-atlas exact-membership pass.
  Sixteen Java formula-reference combinations cover field/PSF and smoothing/
  prior choices, including the lesion class at slot 10 and its output map.
- **P18czm:** the 13-class membership result is handed into base-TOADS rigid
  registration, with MASK excluded and LESION included. The input membership
  bits come from the P18czl Java reference; registration state is compared to
  the already-validated rigid-registration API with those exact reference bits.
  This validates wrapper wiring and class selection; it is not a new Java
  optimizer comparison.
- **P39zzzzv:** the continuing LesionTOADS branch feeds the next thinning,
  refreshed-statistics, energy, and per-voxel relation-aware membership slice.
  The synthetic thinning trace includes 359 accepted and 258 rejected topology
  decisions, and the membership sweep exercises 1,656 direct and 541 prior
  interior voxel/class dispatches.

The run covers the committed code at `d84f203`. All imagery is synthetic. The
full TOADS and LesionTOADS workflows remain incomplete, and this evidence does
not replace the later numerical acceptance and real-data gates.
