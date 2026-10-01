# MaCRUISEplus pre-TGDM to TGDM composition fixture

This fixture contains one project-generated synthetic topology-corrected
initial membership image for the selected TGDM input slot. It uses an
asymmetric folded radial phantom, with foreground membership 200 and
background membership 0. The eight pre-TGDM inputs are loaded from the
committed [`plus_pre_tgdm_composition`](../plus_pre_tgdm_composition/)
fixture, which uses the same 25 x 25 x 25 RAS grid (0.8, 0.8, 1.2 mm).

This is synthetic stage-composition evidence only. It does not test a topology
correction algorithm or imply whole MaCRUISEplus equivalence. The initial
membership is supplied explicitly, as the workflow's topology-corrected seed
must come from an upstream stage not composed here.

Regenerate the image and digest manifest with:

```sh
python tests/data/synthetic/plus_tgdm_composition/generate_fixture.py
```
