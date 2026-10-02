# Synthetic MaCRUISE ACE–GVF–TGDM core composition (P18cybu)

## Ported behavior

The new `run_toads_cruise_core` API and `toads-cruise-core` CLI continue from
second-pass TOADS tissue maps through the selected CRUISE numerical core:
optional MaSeg brainstem/cerebellum removal, GM-derived shared crop,
(6,18) WM-fill topology correction, base ACE, generalized GVF from enhanced
GM, nested TGDM, then inverse crop of inner/central/outer level sets. Base
MaCRUISE leaves the updated multi-atlas CSF ACE option disabled; the
composition explicitly uses that setting. TGDM uses the wrapper's 0.6 WM iso
level, and the selected uncrop restores full-grid output with the source's
background value 20 and exclusive high-bound behavior. Visualization surface
generation is omitted.

## Validation

Focused affected group:

```text
./.venv/bin/pytest -q tests/test_toads_cruise_core.py tests/test_toads_cruise_preparation.py tests/test_toads_tissue_crop.py tests/test_toads_cerebellum.py tests/test_topology.py
```

Result: **26 passed, 0 failed, 0 skipped, 2 expected RAS warnings in 11.44 seconds**.

The committed small asymmetric RAS maps exercise one ACE diffusion iteration,
converged GVF on the compact crop, one iteration for each TGDM profile, the
actual crop-to-core input wiring, native-grid inverse crop, CLI serialization,
and provenance. The core CLI/API outputs are checked for type, RAS geometry,
and source crop padding. An additional default-iteration compact-map check
passed **1 test in 6.18 seconds**; it confirms the actual ACE, GVF, and nested
TGDM default iteration path completes on this synthetic crop. Standalone ACE,
GVF, topology, and TGDM tests provide the component-level references; these
compositions are not a combined Java-wrapper comparison.

This command starts from already-produced second-pass TOADS maps. Historical
SPECTRE registration/transfer, complete default-iteration human-sized
convergence, final surface/thickness products, and the complete MaCRUISE
wrapper remain open.
P18cybu remains Started.
