# P18cg full-suite evidence: populated TOADS grow-then-refresh

- **Implementation commit:** `f7882bb4c54ce72c2bef699153bad15b2a3b780a`
- **Command:** `.venv/bin/pytest -q`
- **Result:** 590 passed, 2 expected RAS reorientation warnings, in 81.94 seconds.
- **Focused check:** `.venv/bin/pytest -q tests/test_toads_growing_pass.py tests/test_toads_growing_queue_initialization.py tests/test_toads_growing_queue_step.py tests/test_toads_initialized_thinning_pass.py tests/test_toads_thinning_queue_run.py tests/test_toads_thinning_queue_step.py tests/test_toads_stopping_iteration.py` — 11 passed in 0.61 seconds.
- **Composition:** the selected Java harness drains the folded 9³ growth queue (448 events), then computes the relation table/map, boundary distances, and closest labels from the grown segmentation. Python's continuing-iteration adapter runs growth even with `doGrowing=false` when the iteration does not stop, and matches each recorded post-growth relation output.
- **Source hash audit:** the correct SHA-256 for the selected VandyCRUISE `MultipleToadSegmentation.java` is `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`. The previously recorded value contained an extra `f`; affected Java fixtures, source-formula harnesses, tests, and two reports were corrected.
- **Limit:** this composes selected growth and relation-refresh methods after a membership update; it does not include membership optimization, Dura removal, registration, or a full TOADS iteration.

See [`tests/data/synthetic/toads_growing_pass/README.md`](../../tests/data/synthetic/toads_growing_pass/README.md) for fixture details and the Java harness generator.
