# In-vivo TGDM stage-isolation attempt — 2026-10-11

## Result

A TGDM-only replay on the saved 161×195×154 RAS crop did not return within a
15-minute cap. It reached the central evolution's signed-distance
reinitialization (`tgdm_nested_arrays` → central `tgdm_evolve` →
`_tgdm_signed_distance_fmm` / FMM march) and was stopped there. The process
ended with exit code 130; a later poll confirmed it had exited. The system
memory snapshot showed 64% free; process-specific RSS was unavailable.

No intermediate NIfTIs/NPYs or stage metrics were produced because
`tgdm_nested_images` returns its intermediate fields only after the full
sequence completes. The saved diagnostic script is under ignored `build/`;
no subject data or partial outputs were committed. Consequently this run does
not identify which TGDM phase changes the ribbon support or component counts.

## Inputs and settings

The diagnostic reuses only the six saved cropped core inputs:

- `wmfill_topology_corrected.nii.gz` as initial membership
- `wmfill_cropped.nii.gz` as white matter
- `gm_cropped.nii.gz` as gray matter
- `ace_enhanced_gm.nii.gz` as ACE gray matter
- `ace_thinned_skeleton.nii.gz` as the skeleton
- `gvf.nii.gz` as the component-last vector field

It also loads the original MaSeg and crops it using the recorded inclusive
bounds `((21,181),(47,241),(88,241))` for fixed support measurements. No
brain-mask argument is supplied, matching the saved core invocation.

Settings match saved provenance: white-matter iso `0.6`; initial/inner/
central/outer iterations `6/4/7/7`; central-shrink iterations `1`; the
remaining parameters use `TGDMBaselineParameters` values. The direct API is
used because the TGDM CLI does not expose this profile and its default
white-matter isovalue differs.

## Reproduction

From the PyMaCRUISE repository root:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python \
  build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/tgdm-stage-isolation-2026-10-11/run_tgdm_stage_isolation.py
```

The failed bounded attempt wrote only that diagnostic script. Its input
inventory and settings are recorded above; the API did not return arrays to
write. This is a real-data runtime diagnostic, not a source-equivalence or
clinical-validation result.

## Next step

Instrument TGDM so each completed phase can persist its field before the next
signed-distance reinitialization, then profile and optimize the central FMM
march against the existing source-selected synthetic comparisons. Retry the
human stage replay under the same memory/runtime cap after the optimized path
passes focused synthetic checks.
