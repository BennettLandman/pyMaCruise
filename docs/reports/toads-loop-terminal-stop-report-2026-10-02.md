# Registered TOADS CLI terminal-condition reporting

The selected registered loop already computed the Java float32 energy ratio
and terminal control decision, but its public result and CLI reported only the
iteration count. The result now exposes the terminal ratio and two independent
flags: whether the energy threshold was reached and whether the configured
iteration limit was reached. Both flags may be true on one terminal iteration,
matching the source order in which either condition selects the stop branch.
The CLI prints these values with the RAS output summary.

## Synthetic cap run

```text
./.venv/bin/pymacruise toads-loop \
  --t1-n4 tests/data/synthetic/toads_registered_first_iteration_cycle/t1.nii.gz \
  --brain-mask tests/data/synthetic/toads_registered_first_iteration_cycle/brain-mask.nii.gz \
  --atlas-descriptor tests/data/synthetic/toads_registered_field_refit/atlas-11.txt \
  --output build/wheel-validation/run/terminal-report-loop-50.nii.gz \
  --registration-levels 1 --registration-iterations 1 \
  --registration-max-iterations 1 --padding 2 \
  --maximum-iterations 50 --relation-limit 1
```

Observed: `iterations=50`, `changed_classification_voxels=189`,
`terminal_energy_ratio=0.547289431`,
`energy_threshold_reached=false`, `iteration_limit_reached=true`. With the
configured energy threshold of `1e-4`, this compact fixture reached the cap
without meeting the energy criterion. Its uint8 output is finite, shape
13×13×13, RAS, and bitwise identical at every voxel to the saved pre-change
packaged-wheel output (zero differing voxels).

The CLI/result change does not alter the segmentation. Focused command
`./.venv/bin/python -m pytest tests/test_toads_loop_cli.py::test_toads_loop_cli_writes_terminal_classification_in_native_ras tests/test_toads_registered_first_iteration_cycle.py tests/test_toads_initial_registration_control.py -q`:
**3 passed, 0 failed, 0 skipped in 49.28s**. Repository-wide full-suite
evidence remains pending.
