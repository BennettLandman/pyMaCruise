# Full-suite evidence — P18cj Dura TGDM evolution

- Date: 2026-09-30
- Tested commit: `03626d4986c20eb9c3fc934d2b09a7ac75f9a84f`
- Dura implementation last changed in: `77e0cfd` (`Compare Dura TGDM against original JIST evolution`)
- Command: `.venv/bin/pytest -q`
- Result: **628 passed**, 2 expected RAS reorientation warnings, in 110.04 seconds.
- Focused validation: 2 passed for the nested full-profile comparison and folded-cortex tolerance gate.

The nested synthetic case matches the selected Java/JIST run's five outer
passes and schedules 150 updates. The angularly folded 25³ cortex case triggers
one source landmine restart; both implementations complete 147 updates. The
final Python output scores Dice **0.997923**, symmetric mean surface distance
**0.005180 mm**, and Hausdorff distance **0.328474 mm**, within the accepted
minimum Dice 0.97, mean distance 0.01 mm, and Hausdorff 1 mm thresholds.

Four voxel signs differ from the Java final field. Intermediate reinitialized
distance fields differ by up to 0.36 voxel in the folded case; these remain
visible measured differences, while the agreed final geometry metrics pass.
The evidence compares the bounded Dura evolution slice and does not establish
full TOADS or full end-to-end MaCRUISE equivalence. The folded fixture and
metrics are described in
`tests/data/synthetic/dura_tgdm_folded_cortex/README.md`.

The two warnings are from established tests that intentionally supply LPS
topology volume and paint mask inputs; both are reoriented to RAS.
