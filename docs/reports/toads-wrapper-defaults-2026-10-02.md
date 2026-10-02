# P18cyah: selected TOADS alignment defaults

The selected `MedicAlgorithmToads08` plugin constructs the TOADS algorithm with
modality default `T1_MPRAGE`, `levels=4`, and `firstIter=50` for initial rigid
alignment. `register_toads_from_atlas_assets` now exposes the same public
defaults while retaining explicit caller overrides for small validation cases.

Source pin: `MedicAlgorithmToads08.java` SHA-256
`d00a7ad609ab0c6381334248b53f4daadc070161af67ef853e2a552dd008d115`.
Implementation commit: `45bd74d885de468d3a5fb6c15303e1729674ea36`.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_atlas_startup.py tests/test_toads_alignment_startup.py tests/test_toads_registration_sampling.py tests/test_toads_registration.py
```

Result on the exact implementation commit: **43 passed, 0 failed, 0 skipped in
2.37 seconds**. The synthetic two-class composition uses explicit low-cost
`levels=2`/`max_iterations=1`, while new assertions inspect the public profile
defaults. This verifies API defaults and override wiring, not execution at the
full four-level/50-iteration setting or full workflow equivalence.
