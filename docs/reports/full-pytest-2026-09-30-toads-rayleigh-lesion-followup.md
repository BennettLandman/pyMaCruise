# Full pytest evidence — TOADS and LesionTOADS fixed-input passes

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest -q` from the repository root |
| Result | **750 passed, 2 warnings** |
| Runtime | 121.12 seconds |
| Warnings | Expected input-orientation warnings: the topology volume and paint mask are reoriented from LPS to RAS in `test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`. |
| Implementation commit | `800e677d52bffb2a1f6c435ddf96b7ac81a2408e` |
| Run date | 2026-09-30 |
| Scope | Complete configured pytest suite for the committed implementation batch. The report/plan documentation commits that followed do not change the tested implementation. |

This run provides full-suite evidence for P18cza, P18czb, and P18czc, plus
P39zzzzn and P39zzzzp. P18cza validates the selected base-TOADS Rayleigh data
energy map; P18czb composes the homogeneous Rayleigh membership pass; P18czc
composes OBJECT/MASK/OUTMASK and preserves the selected Java auxiliary-map
smoothing control flow. P39zzzzn extends the fixed-input LesionTOADS energy
maps to field and special-term branches; P39zzzzp composes the objective-
reporting membership-map pass.

Every new slice has a pinned selected-source Java harness and committed
synthetic RAS NIfTI data. These are selected-method or fixed-input composition
comparisons; they do not establish full base-TOADS, LesionTOADS, or MaCRUISE
workflow equivalence. Fixture provenance is in
`tests/data/synthetic/toads_rayleigh_energy/`,
`tests/data/synthetic/toads_rayleigh_membership/`,
`tests/data/synthetic/toads_outmask_membership/`,
`tests/data/synthetic/lesion_data_energy_maps_special/`, and
`tests/data/synthetic/lesion_objective_membership_pass/`.
