# Full-suite evidence — P18ch and P18ci Dura adapters

- Date: 2026-09-30
- Tested commit: `ce00bbcdfe42c69c8136854b310a3d28d267f149`
- Dura implementation last changed in: `1488bd4` (`Exercise default Dura five-pass wrapper`)
- Command: `.venv/bin/pytest -q`
- Result: **628 passed**, 2 expected RAS reorientation warnings, in 109.34 seconds.
- Focused validation: 2 passed for `test_dura_initial_level_sets_match_original_jist_jar` and
  `test_dura_output_mask_matches_wrapper_threshold_order_and_copies`.

P18ch checks `dura_initial_level_sets`: selected float32 inner/outer membership
level sets and both capped signed-distance volumes match the committed output
from original JIST `DistanceField.class` on the synthetic case. The committed
raw-bit output and generator are
`tests/data/synthetic/dura_initial_levelsets/java_reference.json` and
`scripts/generate_java_dura_initial_levelsets_reference.sh`.

P18ci checks source wrapper comparison order at and around the mask threshold,
level-set translation after mask construction, the resulting T1/WM/GM zeroed
arrays, and the 8-bit mask. It also confirms inputs are not mutated. The
wrapper's full selected stage composition has a separate synthetic RAS
comparison in `test_dura_removal_composes_selected_stages_on_ras_grid`.

The two warnings come from established tests that intentionally supply LPS
topology volume and paint mask inputs; both are reoriented to RAS.
