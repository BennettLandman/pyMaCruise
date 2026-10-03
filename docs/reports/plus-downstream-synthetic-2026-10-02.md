# Synthetic MaCRUISEplus downstream composition

## Executed stages

`tests/test_plus_downstream_integration.py` exercises the selected plus
boundary/tissue/topology/ACE/GVF/TGDM composition followed by post-TGDM
uncropping, surface extraction, and thickness measurement. It uses the
committed 25³ folded-cortex fixture as **precomputed crop-ready atlas and
TOADS products** and an explicit full-grid crop descriptor.

Focused command and result:

```text
./.venv/bin/pytest -q tests/test_plus_downstream_integration.py
1 passed in 14.47s
```

The test checks finite RAS images on the reference geometry, RAS mesh
coordinates in millimeters, finite thickness outputs, expected nested TGDM
products, and that supplied image inputs were not mutated. The composed calls
reuse the existing `run_plus_topology_and_tgdm` and
`run_plus_post_tgdm_stages` APIs.

## Evidence limits and next integration

This is synthetic downstream integration only. It does not execute N4, FLAIR
registration, multi-atlas product generation, crop construction, TOADS, or
LesionTOADS; the fixture supplies products at the crop-ready boundary. It
does not yet provide a restartable CLI/manifest, compare the complete plus
workflow with Java, or validate real imaging data. Next, establish the shared
registered LesionTOADS/atlas-grid input contract and turn this downstream
composition into checked, restartable workflow stages.
