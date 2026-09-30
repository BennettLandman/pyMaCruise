# Full pytest evidence: LesionTOADS final classification

- **Implementation commit:** `2b7515c415a6f4f65371c279d611ca4c50f4c5a3`
- **Command:** `.venv/bin/pytest -q`
- **Result:** 598 passed, 2 expected RAS-reorientation warnings, 87.88 seconds.
- **Focused selection:** `tests/test_lesion_toads_final_classification.py` — 4 passed.

`lesion_toads_final_classification` composes the selected final label transfer
with `updateLesionClassifcation`. The Java selected-method excerpt and Python
match the mapped class labels, changed-label count, full 0/1/2 lesion subtype
map, and lesion voxel count on the committed 3×2×2 synthetic RAS volumes.
Inputs exercise strict membership comparisons, equal lesion/black-hole
memberships, non-lesion labels, and `hasLesion=false`. The membership and
classification input NIfTIs are stored as `.nii.gz`; their axis codes and
1-mm voxel sizes are checked.

This validates output label transfer and subtype assignment for the supplied
memberships only. It does not validate upstream LesionTOADS membership or
topology evolution, plus-pipeline outputs, or clinical/equivalence metrics.
