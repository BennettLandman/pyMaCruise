# LesionTOADS registered loop driver

- Focused command: `./.venv/bin/pytest -q tests/test_lesion_registered_repeated_iteration.py`
- Focused result: **1 passed, 0 failed, 0 skipped in 9.45s**.
- Tested implementation commit: recorded after this implementation slice is committed.

Added `run_lesion_toads_descriptor_loop` to carry the selected 13-class branch
through repeated thinning, statistics, membership, energy-ratio, registration,
shape resampling, and source stop control. The asymmetric 9³ RAS integration
compares its second and third updates and final iteration-limit classification
with the same stages composed directly. Shape maps use the registration's
updated support bounds at each pass.

The caller supplies initial registration and wrapper state. The selected
constructor distance maps and correction fields are fixed across this bounded
loop; dynamic distance refresh, complete initialization/preprocessing,
unconstrained convergence, complete LesionTOADS, and real MRI remain outside
this adapter. The synthetic composition is not a composite Java equivalence
result. Full-suite evidence is pending.
