# MaCRUISEplus entry-point option contract

This note maps the selected MATLAB `maCRUISEpluspipeline.m` option overrides
to a Python settings object. It is a configuration record only; it does not
claim that the ACE, topology, CRUISE, or LesionTOADS computations have been
ported by the settings object.

## Selected option values

| MATLAB option | Selected plus value | Python field | Use in the selected pipeline |
|---|---:|---|---|
| `N3orN4` | `N4` | `bias_correction` on `WorkflowConfig` | Correct T1 and FLAIR before the lesion preparation branch. |
| `LesionToadsOutput` | `dura removal inputs` | `lesion_toads_output` | Request the LesionTOADS products consumed by dura removal and the later plus path. |
| `LesionToadsCorrect` | `false` | `lesion_toads_correct` | Disable LesionTOADS's own field correction. |
| `ErodeMethod` / `OverLookLabels` | `3D` / `120,121,202,203,156,157` | `boundary_erosion_method` / `boundary_overlook_labels` | Select 3D boundary erosion and labels passed to boundary construction. `run_plus_boundary_stage` constructs the boundary, then removes voxels where aligned TOADS WM is greater than 0.01; it returns raw and cleaned products. Synthetic routing evidence: [P40g full-suite report](reports/full-pytest-2026-09-29-plus-boundary-routing.md). |
| `update_Enhance` | `0` | `update_enhance` | Disable enhancement update. |
| `WMIsoLevel` | `0.5` | `wm_iso_level` | `plus_tgdm_parameters` maps this to `TGDMBaselineParameters.white_matter_iso` while preserving the remaining baseline values. See the [P40d parameter-routing report](reports/full-pytest-2026-09-29-plus-wm-isolevel-routing.md). |
| `GVF_Lambda` | `0.2` | `gvf_lambda` | `run_plus_gvf_stage` forwards this value to `generalized_gvf_image` with the source wrapper's 200-iteration maximum and no vector normalization. See the [P40c parameter-routing report](reports/full-pytest-2026-09-29-plus-gvf-routing.md). |
| `IfGenerateIntermediate` | `true` | `generate_intermediate` | Request intermediate output files in the selected source workflow. |
| `NewCSFfile` | output of boundary cleaning | not a scalar option | Supplies the updated CSF/boundary input to ACE. |
| `UseUpdatedCSFinACE` | `3` | `use_updated_csf_in_ace` | `run_plus_ace_stage` passes placement 3 to the existing ACE API. Merge positive updated-CSF voxels after the speed map's second alpha-1.8 Gaussian; the selected MaACE wrapper writes speed 0.1 there. See [the P40b Java comparison](reports/full-pytest-2026-09-29-plus-ace-after-gaussian-csf-merge.md) and [the P40e wiring report](reports/full-pytest-2026-09-29-plus-maace-stage.md). |
| `MaWM` | cropped pre-update WM file | not a scalar option | Input used for topology refinement. |
| `refineTopoCorrectByMa` | `1` | `refine_topology_by_maseg` | `run_plus_topology_stage` requires and forwards MaWM into the existing atlas-refinement branch; with the flag off, it preserves base topology output. See the [P40f routing report](reports/full-pytest-2026-09-29-plus-topology-routing.md). |

The original plus entry point applies these values after image preparation,
the boundary generation/cleaning step, and the GM/WM update. The Python
`PlusPipelineOptions` captures only scalar/boolean option values; actual file
products remain explicit stage inputs so their raw, corrected, dura-processed,
and cropped provenance can be represented by the eventual runner.

For the GM/WM update, the selected wrapper calls the Java
`CombineMAandTOADS` class with `ErodeMethod=3D` and `GMconfident=0.001`.
`run_plus_tissue_update_stage` routes the six aligned cropped atlas/TOADS
products to the existing source-mapped Python combiner and returns its GM, WM,
and mask images. Synthetic routing evidence is recorded in the
[P40h report](reports/full-pytest-2026-09-29-plus-tissue-update.md).

The selected surface wrapper supplies the cropped corrected level set,
updated WM/GM, MaACE-enhanced GM and thinned skeleton, and the plus GVF field to
nested TGDM. `run_plus_tgdm_stage` maps WMIsoLevel into the existing nested
image adapter and returns level sets on that cropped RAS grid; see the
[P40i report](reports/full-pytest-2026-09-29-plus-tgdm-routing.md).
`uncrop_plus_tgdm_outputs` restores its six scalar volume products using the
same GM-derived crop bounds and selected JIST exclusive-maximum convention;
this does not convert mesh vertices. See the
[P40j report](reports/full-pytest-2026-09-29-plus-tgdm-uncrop.md).

## Image provenance and explicit input mask

The selected MATLAB branch sends T1 and FLAIR through N4 before stripping,
segmentation, and dura processing. It later passes the original T1 path to
the surface pipeline alongside prepared GM and updated/cropped tissue maps.
PyMaCRUISE requires the caller's brain mask on the T1 grid. That input
supplies the mask needed by Python preparation in place of the historical
SPECTRE/BET skull-strip step; no BET or FreeSurfer dependency is introduced.
This does not imply that all original/corrected image handoffs are equivalent;
the stage runner must preserve and test them explicitly.

## Source files

- `pipeline/maCRUISEpluspipeline.m` selects N4, the LesionTOADS outputs and
  correction mode, boundary behavior, and surface-stage overrides.
- `pipeline/run_LesionToads_pipeline.m` shows the N4 → mask/strip → EM → dura
  → lesion-TOADS ordering and its image handoffs.
- `pipeline/run_LesionToads_jist.m` records the selected LesionTOADS atlas,
  output and convergence arguments.

The code mapping is in `src/pymacruise/plus_pipeline.py`. The base workflow
does not consume `PlusPipelineOptions`, so these plus-specific overrides do
not alter its settings.
