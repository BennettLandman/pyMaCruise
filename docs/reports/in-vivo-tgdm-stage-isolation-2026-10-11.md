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

## Compiled-FMM rerun and held-fixed Java source comparison

The same TGDM-only replay completed after the source-ordered Numba FMM
optimization. On the saved 161×195×154 human crop it took **315.5 s** and
peaked at **3.418 GiB RSS** under the 12 GiB watchdog. The terminal inner and
outer surfaces support **81.963%** of input MaSeg GM and **83.253%** of input
MaSeg WM, so this run still fails the agreed 97% support criterion. It saved
all seven phase fields and support/component metrics under the ignored
`build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/tgdm-stage-isolation-2026-10-11/`
directory. The first recorded TGDM membership already has 209 6-connected
components; subsequent fields keep 209, so TGDM does not initiate the
observed 6-connected component count.

A held-fixed Java replay using the MaCRUISE-specific `MaGenericTGDM` source
and the selected crop-wrapper settings matches all five surface sign masks
exactly (Dice 1.0). Its final outer field has Dice **0.9999792**, symmetric
mean boundary-voxel distance **0.000333 mm**, and Hausdorff **0.5764 mm**.
These are within the agreed limits, but are 6-neighbor voxel-center
approximations, not the formal mesh-to-mesh metric. The output fields show
small numerical differences; see the ignored comparison JSON for per-stage
MAE and maxima.

A separate VandyCRUISE `GenericTGDM` source-profile replay was run to honor the
requested Vandy source audit. It matches the same fields through central
evolution, then diverges: central-shrink Dice **0.9896**, final inner Dice
**0.9422**, and outer Dice **0.9968**. Approximate boundary-shell mean/maximum
distances are 0.094/1.414 mm, 0.329/7.280 mm, and 0.037/4.359 mm, respectively.
These shell metrics are diagnostic approximations, not mesh distances. Source
inspection found that Vandy `GenericTGDM.updatePoint` leaves persistent
advection state untouched when `advectionForces` is null, whereas MaCRUISE's
`MaGenericTGDM` explicitly clears it. This difference begins after central
evolution uses a non-null GVF. The crop wrapper selects the MaCRUISE-specific
class, so the Vandy generic-class mismatch is documented as a source-variant
risk and is not used to alter the Python implementation. Both Java replays
completed below their memory caps (2.204 GiB for MaGeneric; 1.627 GiB for
Vandy Generic).

The saved per-stage Python support report is diagnostic: it pairs each inner
candidate with the final outer surface and vice versa. The formal terminal
support remains the 81.963% GM / 83.253% WM result. This is an in-vivo,
held-fixed component comparison, not full-pipeline source equivalence or
clinical validation.

## Next step

Continue investigation at the terminal shell/label support handoff and compare
the MaGeneric phase fields with the saved topology-corrected inputs. Do not
change source parameters or relax support thresholds without evidence from the
selected MATLAB/Java wrapper.
