# Full pytest evidence — TOADS Rayleigh and LesionTOADS energy batch

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest -q` from the repository root |
| Result | **736 passed, 2 warnings** |
| Runtime | 168.07 seconds |
| Warnings | Expected input-orientation warnings: the topology volume and paint mask are reoriented from LPS to RAS in `test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`. |
| Implementation commit | `a203b7ec912046b07017de34b63cb595a5f887af` |
| Run date | 2026-09-30 |
| Scope | Complete configured pytest suite at the implementation commit above. |

This run provides the requested full-suite evidence for P18cy, P18cz, P39x,
and P39zzzzk–P39zzzzm. P18cy's selected Java dispatch harness matches the
interleaved OBJECT/MASK/optimized tissue memberships and objective. P18cz's
Rayleigh method excerpt matches all 17 IEEE edge cases bitwise. P39x uses the
corrected selected-source Sulcal-CSF control flow: unfielded execution adds
the T1 special term followed by generic weighted channels, while fielded
execution makes the T1 term exclusive. P39zzzzk, P39zzzzl, and P39zzzzm
independently validate the fixed-input LesionTOADS membership-map pass,
objective reduction, and raw volume energy maps.

These are component comparisons against pinned Java excerpts/harnesses and
synthetic RAS NIfTI fixtures. They do not establish full base-TOADS or
LesionTOADS workflow equivalence. The corresponding fixture provenance is in
`tests/data/synthetic/toads_optimized_dispatch/`,
`tests/data/synthetic/toads_rayleigh_distance/`,
`tests/data/synthetic/lesion_toads_sulcal_csf/`,
`tests/data/synthetic/lesion_membership_map_pass/`,
`tests/data/synthetic/lesion_objective_reduction/`, and
`tests/data/synthetic/lesion_data_energy_maps/`.
