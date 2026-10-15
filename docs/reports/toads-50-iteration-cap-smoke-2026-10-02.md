# Registered TOADS 50-iteration CLI smoke check (P18cyby)

## Run

The built wheel from [P18cybn](toads-loop-wheel-2026-10-02.md) ran
`toads-loop` with the committed 13³ synthetic T1, brain mask, and 11-class
atlas. Registration was held to one level and one iteration to isolate the
outer registered TOADS loop; the loop limit was set to the source default of
50, with relation limit 1 and padding 2.

```text
PYTHONPATH=build/wheel-validation/site build/wheel-validation/site/bin/pymacruise toads-loop \
  --t1-n4 tests/data/synthetic/toads_registered_first_iteration_cycle/t1.nii.gz \
  --brain-mask tests/data/synthetic/toads_registered_first_iteration_cycle/brain-mask.nii.gz \
  --atlas-descriptor tests/data/synthetic/toads_registered_field_refit/atlas-11.txt \
  --output build/wheel-validation/run/default-loop.nii.gz \
  --registration-levels 1 --registration-iterations 1 \
  --registration-max-iterations 1 --padding 2 \
  --maximum-iterations 50 --relation-limit 1
```

Initial packaged-wheel result: `iterations=50`,
`changed_classification_voxels=189`. A source-checkout rerun after the CLI
gained stop telemetry reports terminal ratio `0.547289431`,
`energy_threshold_reached=false`, and `iteration_limit_reached=true` for the
same inputs and settings. See the later [terminal-condition report](toads-loop-terminal-stop-report-2026-10-02.md).

## Output check

The wheel-produced classification is uint8, shape 13×13×13, RAS, finite, and
contains all 11 expected descriptor labels (1–11). It completed the full
configured loop cap on the compact case and is bitwise identical to the
pre-change packaged-wheel output (zero differing voxels). The later
source-checkout report confirms this is an iteration-cap stop, not natural
convergence. This is not Java-wrapper equivalence; the follow-up telemetry
change and its focused/full-suite validation are recorded in the linked
terminal-condition report.

The multi-minute run on a 13³ phantom is a performance risk to quantify before
human-sized execution. The report does not extrapolate its duration linearly to
human images.
