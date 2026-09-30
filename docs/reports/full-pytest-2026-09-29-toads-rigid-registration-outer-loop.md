# P18cc full-suite evidence: repeated TOADS rigid registration updates

- **Implementation commit:** `2df56c29bfb6f880e893be28c8333098591c175b`
- **Command:** `.venv/bin/pytest -q`
- **Result:** 584 passed, 2 expected RAS reorientation warnings, in 86.52 seconds.
- **Focused check:** `.venv/bin/pytest -q tests/test_toads_rigid_registration_update.py` — 2 passed in 0.38 seconds.
- **Reference:** the committed selected-source Java formula harness runs the rigid registration outer loop on a deterministic synthetic 4×4×4 membership volume and 5×5×5 three-class atlas. Two updates improve the objective; the last rejects 11 trials and stops. Python matches each iteration's transform, energy, lambda, relative change, acceptance state, and rejection count.
- **Limit:** this checks the selected optimizer loop and coefficient update on one small synthetic case. It does not validate atlas-to-subject alignment, pyramid setup, deformable TOADS, segmentation, or clinical behavior.

The fixture and regeneration details are in [`tests/data/synthetic/toads_rigid_registration_update/README.md`](../../tests/data/synthetic/toads_rigid_registration_update/README.md).
