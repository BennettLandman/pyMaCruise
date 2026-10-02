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

CLI result: `iterations=50`, `changed_classification_voxels=189`.

## Output check

The wheel-produced classification is uint8, shape 13×13×13, RAS, finite, and
contains all 11 expected descriptor labels (1–11). It completed the full
configured loop cap on the compact case. The CLI result does not expose a
separate stop reason, so reaching iteration 50 is not evidence that the
convergence threshold triggered. This is an exploratory wheel smoke check,
not a pytest or Java-wrapper equivalence result; it is recorded as Started in
`plan.md` pending a repeatable test and full-suite checkpoint.

The multi-minute run on a 13³ phantom is a performance risk to quantify before
human-sized execution. The report does not extrapolate its duration linearly to
human images.
